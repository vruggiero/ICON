!
!+ Fortran 90 Interfaces to the  ECMWF EMOS GRIB library or MPIfM clone
!
MODULE mo_emos_grib1
!
! Description:
!   Fortran 90 Interfaces to the GRIBEX routine
!              (ECMWF EMOS library or MPI clone)
!
!   A data type 't_grib1' is defined, to hold the content of a GRIB-record
!   (both encoded and decoded data) and some information of the associated
!   file (name, position, record number). The decoded sections are stored
!   in components of derived types, so that specific entries can be
!   specified by name (and not by position within an integer
!   array). Depending on the grid representation, different data types are
!   provided for section 2 (Grid Definition Block).
!
!   Wrapper routines are provided for the ECMWF EMOS library
!   routines. These wrapper routines only take a variable of type
!   't_grib1' as their argument. All input/output from/to the library
!   routines is stored within this variable. Array components for the
!   encoded and decoded data are allocated with the correct size.
!
!   A GRIB2 record handle is provided for the case that GRIB decoding
!   is performed by the GRIB2 API. If this handle is valid the routines
!   and derived types of this module are bypassed.
!
!   Contents of ECMWF EMOS library:
!   (* denotes wrapper routine defined in this module)
!   !------------+-------------------+-----------------------------------------
!   !   Fortran  | C interface       | Description
!   !------------+-------------------+-----------------------------------------
!   ! * GRIBEX   | gribExDP          | Encode/Decode GRIB record
!   !   GRPRS0   | gribPrintSec0     | Print section 0 (IS)
!   !   GRPRS1   | gribPrintSec1     | Print section 1 (PDS)
!   !   GRPRS2   | gribPrintSec2DP   | Print section 2 (GDS)
!   !   GRPRS3   | gribPrintSec3DP   | Print section 3 (BMS)
!   !   GRPRS4   | gribPrintSec4DP   | Print section 4 (BDS)
!   !   GRPRS4W  | gribPrintSec4Wave | Print the wave coordinate information
!   !   GRSDBG   | gribSetDebug      | Set debug printing
!   !   GRSRND   | gribSetRound      | Set rounding
!   !   GRSREF   | gribSetRefDP      | Set reference value
!   !------------+-------------------+-----------------------------------------
!   ! * PBOPEN   | gribOpen          | Open a GRIB file
!   ! * PBCLOSE  | gribClose         | Close a GRIB file
!   !   PBREAD   | gribRead          | Read a block of bytes
!   ! * PBWRITE  | gribWrite         | Write a block of bytes
!   ! * PBSEEK   |                   | Seeks to a specified location
!   !   PBTELL   | gribGetPos        | Tells current byte offset
!   !   PBFLUSH  |                   | Flushes file
!   !   PBSIZE   | gribGetSize       | Size in bytes of the next GRIB product
!   ! * PBGRIB   |                   | Read next GRIB product
!   !------------+-------------------+-----------------------------------------
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
! V1_4         2009/03/26 Harald Anlauf
!  fixes, changes for NEC SX
! V1_5         2009/05/25 Andreas Rhodin
!  work around bug in MPIfM GRIB library (run type/exp)
! V1_6         2009/06/10 Harald Anlauf
!  gribex_90_2: optional argument 'scanmode'
! V1_7         2009/08/24 Andreas Rhodin
!  diff_grib: print octets present in only one of the GRIBs to compare
! V1_8         2009/12/09 Andreas Rhodin
!  cope with different size of local extension 253 by DWD and MPI GRIBlib
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  increase size of isec2 from 1024 to 4096 for operational IFS Gaussian grid
! V1_13        2011/11/01 Andreas Rhodin
!  optimisations, changes for GRIB2 API
! V1_15        2011/12/06 Andreas Rhodin
!  handle the case that all data has the same value and is not encoded on input
! V1_19        2012-04-16 Harald Anlauf
!  pbopen: print 'mode' when opening file fails
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON; read/write proper ICON grid metadata.
!  WMO code, Temporary hack for different size of GRIB1/Sect.1 for MeteoSwiss.
!  FTRACE instrumentation
! V1_27        2013-11-08 Daniel Leuenberger
!  Fixes for GRIB1 local extensions for centre=Meteoswiss,COSMO
! V1_42        2015-06-08 Harald Anlauf
!  ICON local patch
! V1_46        2016-02-05 Harald Anlauf
!  Handle DWD local definiton 252, GRIB2 encoding of mean and spread
! V1_47        2016-06-06 Andreas Rhodin
!  work around wrong value returned by PBSIZE
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2003-2008  original code
! Harald Anlauf   DWD  2007-2008  bug fixes, optimisations, extensions
!------------------------------------------------------------------------------
! In order to get interp_strato working in ICON, the _DACE_ in the following line
! should be replaced by __DACE__:
#if defined (__ICON__) && !defined (_DACE_)
#define NO_GRIBEX
#endif
!------------------------------------------------------------------------
! On IBM the ichar intrinsic returns negative numbers.
! (The treatment of character codes > 127 is not defined by the standard)
! Thus we provide this replacement.
!------------------------------------------------------------------------
#if defined(__ibm__)
#define ICHAR jchar
#endif
!------------------------------------------------------------------------------
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && defined(NOMPI)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!------------------------------------------------------------------------------
#include "tr15581.incf"
  !=============
  ! Modules used
  !=============
  USE mo_kind,       ONLY: sp, dp, wp, i8  ! kind parameters
  USE mo_exception,  ONLY: finish, message ! abort on error condition
  USE mo_mpi_dace,   ONLY: p_bcast,       &! generic broadcast routine
                           p_bcast_ptr,   &! pointer broadcast routine
                           dace            ! default communicator
  USE mo_ecmwf_grib, ONLY: t_s1_ecmwf,    &! ecmwf local extension data type
                           decode_ecmwf,  &! decode the ECMWF local extensions
                           encode_ecmwf,  &! encode the ECMWF local extensions
                           grprs1ec        ! print ECMWF local use part of sec1
  USE mo_endian,     ONLY: little          ! returns .true. for little endian
  USE mo_wmo_tables
  IMPLICIT NONE

  !================
  ! Public entities
  !================
  PRIVATE
  !----------------------
  ! Data type definitions
  !----------------------
  PUBLIC :: t_grib1         ! GRIB record data type
  PUBLIC :: t_sec0          ! component: Section0
  PUBLIC :: t_sec1          !            Section1, product definition.
  PUBLIC :: t_s1_dwd        !            Section1, DWD specific section.
  PUBLIC :: t_rsec2         !            Section2, real values in grid def.
  PUBLIC :: t_s2_gauss      !            Section2, Gaussian grid.
  PUBLIC :: t_s2_latlon     !            Section2, lat/lon grid.
  PUBLIC :: t_s2_tri        !            Section2, triangular grid.
  PUBLIC :: t_s2_sph        !            Section2, spherical harmonics.
  PUBLIC :: t_sec3          !            Section3, bitmap (integer field).
  PUBLIC :: t_rsec3         !            Section3, bitmap (real field).
  PUBLIC :: t_sec4          !            Section4, binary data
  !------------------------------------
  ! Fortran 90 interface to a subset of
  ! the EMOS GRIB and BPIO library
  !------------------------------------
  PUBLIC :: pbopen          ! Open a GRIB file
  PUBLIC :: pbclose         ! Close a GRIB file
  PUBLIC :: pbgrib          ! Read next GRIB product
  PUBLIC :: pbwrite         ! Write a block of bytes
  PUBLIC :: pbseek          ! Seeks to a specified location
  PUBLIC :: gribex          ! Encode/Decode GRIB record
  PUBLIC :: pbsetraw        ! Set DWD record headers (or raw mode)
  PUBLIC :: grprs1lu        ! Print local use part of section 1
#ifdef USE_PBSEEK64
  PUBLIC :: pbtell          ! Return current file pointer
#endif
  !-------------------------------------------
  ! additional operations on data type t_grib1
  !-------------------------------------------
  PUBLIC :: reallocate_data   ! (re)allocate  t_grib1% rsec4
  PUBLIC :: reallocate_buffer ! (re)allocate  t_grib1% kgrib
  PUBLIC :: assignment (=)    ! copy GRIB record
  PUBLIC :: p_bcast           ! generic broadcast routine
  PUBLIC :: bufr_size_enc     ! determine size of encoded BUFR
  PUBLIC :: bufr_size_dec     ! determine size of decoded BUFR
  PUBLIC :: diff_grib         ! compares 2 grib-records
  PUBLIC :: local_dwd         ! DWD   local extension present ?
  PUBLIC :: local_ecmwf       ! ECMWF local extension present ?
  PUBLIC :: get_octets        ! decode octets
  PUBLIC :: put_octets        ! encode octets
  PUBLIC :: print_isec4       ! print integer data of section 4
  PUBLIC :: destruct          ! deallocate components
  !------------------------
  ! extension for GRIB2 API
  !------------------------
  PUBLIC :: INVALID_HANDLE    ! invalid value for GRIB-API handle

  !============================
  ! GRIB 1 data type definition
  !============================
  !----------
  ! Section 0
  !----------
  TYPE t_sec0
    SEQUENCE
    INTEGER :: n_octets  ! 1 - Number of octets in the GRIB message
    INTEGER :: edition   ! 2 - GRIB edition number
  END TYPE t_sec0

  !---------------------------------------
  ! Section 1 (product definition section)
  !---------------------------------------
  TYPE t_sec1
    SEQUENCE
    INTEGER :: table           !  1 | Version number of code table 2 (ECMWF local code table 2).
    INTEGER :: center          !  2 | Identification of centre       (WMO code table 0).
    INTEGER :: process         !  3 | Generating process identification number.
    INTEGER :: grid            !  4 | Grid definition.
    INTEGER :: present_2_3     !  5 | Bit field. 128: Section 2 is included; 64: Section 3 is included.
    INTEGER :: code            !  6 | Parameter indicator (WMO code table 2).
    INTEGER :: level_type      !  7 | Type of level indicator (see WMO code table 3) or satellite identifier.
    INTEGER :: level_st        !  8 | Height, pressure, etc of level (WMO code table 3). Single level or top of layer.
    INTEGER :: level_b         !  9 | Height, pressure, etc of level (see WMO code table 3). Bottom of layer
    INTEGER :: year            ! 10 | Year of century (YY) of Reference time of data. (eg 1 to 100 for year 1901 to 2000).
    INTEGER :: month           ! 11 | Month           (MM) of Reference time of data.
    INTEGER :: day             ! 12 | Day             (DD) of Reference time of data.
    INTEGER :: hour            ! 13 | Hour            (HH) of Reference time of data.
    INTEGER :: minute          ! 14 | Minute          (MM) of Reference time of data.
    INTEGER :: time_unit       ! 15 | Time unit indicator (WMO code table 4).
    INTEGER :: p1              ! 16 | Time period (number of time units) 0 for analyses.
    INTEGER :: p2              ! 17 | Time period. Or time interval between analyses, or forecasts undergoing.
    INTEGER :: time_range      ! 18 | Time range indicator (WMO code table 5).
    INTEGER :: n_average       ! 19 | Number of products included in an average.
    INTEGER :: n_missing       ! 20 | Number of products missing from an average.
    INTEGER :: century         ! 21 | Century         (CC) of Reference time of data. (eg 20 for years 1901 to 2000).
    INTEGER :: sub_center      ! 22 | Sub-centre identifier.
    INTEGER :: factor          ! 23 | Decimal scale factor.
    INTEGER :: local_flag      ! 24 | indicate local use in Section 1. 0: No local use, 1: Local use
    INTEGER :: reserved(25:36) !    | Reserved for WMO reserved fields. Set to 0.
    INTEGER :: local_ident     ! 37 | ECMWF local GRIB use definition identifier. (ECMWF local GRIB usage definitions).
    INTEGER :: local_use   (38:1024)! used for ECMWF local extensions. 192 to 255 free for use by Member States.
  END TYPE t_sec1

  !---------------------------------------
  ! Section 1 (product definition section)
  ! DWD specific part
  !---------------------------------------
  TYPE t_s1_dwd                ! octets|
    SEQUENCE
    INTEGER :: local_ident     ! 41    | local GRIB use identifier.       (254)
    INTEGER :: day_number      ! 42    | not used any more                (255)
    INTEGER :: record_number   ! 43-45 | not used any more                (200)
    INTEGER :: decoding        ! 46    | not used any more                (255)
    INTEGER :: element_no      ! 47    | not used any more                (255)
    INTEGER :: year            ! 48    | year  (without century)
    INTEGER :: month           ! 49    | month (creation date..)
    INTEGER :: day             ! 50    | day
    INTEGER :: hour            ! 51    | hour
    INTEGER :: minute          ! 52    | minute
    INTEGER :: exp             ! 53-54 ! experiment no
    INTEGER :: run_type        ! 53-54 ! run_type (* 2**14 + exp)
    INTEGER :: user_id         ! 55    ! User id, specified by table
    INTEGER :: experiment_id   ! 56-57 ! Experiment identifier
    INTEGER :: ensemble_id     ! 58-59 ! Ensemble identification by table
    INTEGER :: ensemble_size   ! 60-61 ! Number of ensemble members
    INTEGER :: ensemble_no     ! 62-63 ! Actual number of ensemble member
    INTEGER :: major_version   ! 64    ! Model major version number
    INTEGER :: minor_version   ! 65    ! Model minor version number
  END TYPE t_s1_dwd

  !------------------------------------
  ! Section 2 (grid definition section)
  ! Triangular grids
  !------------------------------------
  TYPE t_s2_tri
    SEQUENCE
    INTEGER :: repr        !  1 | 192  Data representation type (table 6).
    INTEGER :: ni2         !  2 | Number of factor 2 in factorisation of NI.
    INTEGER :: ni3         !  3 | Number of factor 3 in factorisation of NI.
    INTEGER :: nd          !  4 | Number of diamonds.
    INTEGER :: ni          !  5 | Number of triangular subdivisions.
    INTEGER :: orient      !  6 | Flag for orientation of diamonds.
    INTEGER :: lat_pole    !  7 | Latitude of pole point.
    INTEGER :: lon_pole    !  8 | Longitude of pole point.
    INTEGER :: lon_dia1    !  9 | Longitude of the first diamond.
    INTEGER :: scan_mode   ! 10 | Flag for storage sequence.
    INTEGER :: reserved_11 ! 11 |
    INTEGER :: nvcp        ! 12 ! Number of vertical coordinate parameters.
!   INTEGER :: reserved (13: 22)! reserved.
    !-----------------------------------------------------------
    INTEGER :: reserved (13: 15)! reserved.
    INTEGER :: npts        ! 16 ! numberOfDataPoints
    INTEGER :: grid_num    ! 17 ! numberOfGridUsed
    INTEGER :: grid_ref    ! 18 ! numberOfGridInReference
    CHARACTER :: uuid(16) !19-22! ICON unique horizontal grid ID
    !-----------------------------------------------------------
  END TYPE t_s2_tri

  !-----------------------------------------------------------------
  ! Section 2 (grid definition section)
  ! Latitude/longitude or equidistant cylindrical or plate carree grids
  !-----------------------------------------------------------------
  TYPE t_s2_latlon
    SEQUENCE
    INTEGER :: repr        !  1 | 192  Data representation type (WMO code table 6).
    INTEGER :: ni          !  2 | Number of points along a parallel.
    INTEGER :: nj          !  3 | Number of points along a meridian.
    INTEGER :: lat_first   !  4 | Latitude of the first grid point.
    INTEGER :: lon_first   !  5 | Longitude of the first grid point.
    INTEGER :: increments  !  6 | 128: Direction increments given, 0: not given.
    INTEGER :: lat_last    !  7 | Latitude of the last grid point.
    INTEGER :: lon_last    !  8 | Longitude of the last grid point.
    INTEGER :: di          !  9 | i direction increment.
    INTEGER :: dj          ! 10 | j direction increment.
    INTEGER :: scan_mode   ! 11 | Scanning mode flags (WMO code table 8).
    INTEGER :: nvcp        ! 12 | Number of vertical coordinate parameters.
    INTEGER :: lat_rot     ! 13 | Latitude of the southern pole of rotation.
    INTEGER :: lon_rot     ! 14 | Longitude of the southern pole of rotation.
    INTEGER :: lat_strech  ! 15 | Latitude of the pole of stretching.
    INTEGER :: lon_strech  ! 16 | Longitude of the pole of stretching.
    INTEGER :: reduced     ! 17 | 0: Regular grid, 1: Quasi-regular (reduced) grid.
    INTEGER :: earth       ! 18 | 0: spherical r=6367.47km, 64: oblate spheroidal (IAU1965)
    INTEGER :: components  ! 19 | 0: u and v components in E/N, 8: components relative to grid
    INTEGER :: reserved (20:22) ! reserved.
    !                    23-nn  | Specification of quasi-regular (reduced) grid points.
  END TYPE t_s2_latlon

  !------------------------------------
  ! Section 2 (grid definition section)
  ! Gaussian grids
  !------------------------------------
  TYPE t_s2_gauss
    SEQUENCE
    INTEGER :: repr        !  1 | Data representation type (WMO code table 6).
    INTEGER :: ni          !  2 | Number of points along a parallel.
    INTEGER :: nj          !  3 | Number of points along a meridian.
    INTEGER :: lat_first   !  4 | Latitude of the first grid point.
    INTEGER :: lon_first   !  5 | Longitude of the first grid point.
    INTEGER :: increments  !  6 | 128: Direction increments given, 0: not given.
    INTEGER :: lat_last    !  7 | Latitude of the last grid point.
    INTEGER :: lon_last    !  8 | Longitude of the last grid point.
    INTEGER :: di          !  9 | i direction increment.
    INTEGER :: nglh        ! 10 | Number of parallels between a pole and the Equator.
    INTEGER :: scan_mode   ! 11 | Scanning mode flags (WMO code table 8).
    INTEGER :: nvcp        ! 12 | Number of vertical coordinate parameters.
    INTEGER :: lat_rot     ! 13 | Latitude of the southern pole of rotation.
    INTEGER :: lon_rot     ! 14 | Longitude of the southern pole of rotation.
    INTEGER :: lat_strech  ! 15 | Latitude of the pole of stretching.
    INTEGER :: lon_strech  ! 16 | Longitude of the pole of stretching.
    INTEGER :: reduced     ! 17 | 0: Regular grid, 1: Quasi-regular (reduced) grid.
    INTEGER :: earth       ! 18 | 0: spherical r=6367.47km, 64: oblate spheroidal (IAU1965)
    INTEGER :: components  ! 19 | 0: u and v components in E/N, 8: components relative to grid
    INTEGER :: reserved (20:22) ! reserved.
    !                    23-nn  | Specification of quasi-regular (reduced) grid points.
  END TYPE t_s2_gauss

  !------------------------------------
  ! Section 2 (grid definition section)
  ! Spherical harmonic coefficients
  !------------------------------------
  TYPE t_s2_sph
    SEQUENCE
    INTEGER :: repr        !  1 | Data representation type (WMO code table 6).
    INTEGER :: j           !  2 | J pentagonal resolution parameter.
    INTEGER :: k           !  3 | K pentagonal resolution parameter.
    INTEGER :: m           !  4 | M pentagonal resolution parameter.
    INTEGER :: repr_type   !  5 | 1 Representation type (WMO code table 9)
    INTEGER :: repr_mode   !  6 | 1: normal packing, 2: complex packing.
    INTEGER :: reserved1 ( 7:11)! Reserved. Set to 0.
    INTEGER :: nvcp        ! 12 | Number of vertical coordinate parameters.
    INTEGER :: lat_rot     ! 13 | Latitude of the southern pole of rotation.
    INTEGER :: lon_rot     ! 14 | Longitude of the southern pole of rotation.
    INTEGER :: lat_strech  ! 15 | Latitude of the pole of stretching.
    INTEGER :: lon_strech  ! 16 | Longitude of the pole of stretching.
    INTEGER :: reserved2 (17:22)! Reserved. Set to 0.
  END TYPE t_s2_sph

  !------------------------------------
  ! Section 2 (grid definition section)
  ! real numbers
  !------------------------------------
  TYPE t_rsec2
    SEQUENCE
    REAL(dp) :: rot_angle  !  1 | Angle of rotation.
    REAL(dp) :: str_factor !  2 | Stretching factor.
    REAL(dp) :: reserved (3: 10)! Reserved. Set to 0.
    REAL(dp) :: vcp      (1:502)! (11:512) Vertical coordinate parameters.
  END TYPE t_rsec2

  !--------------------------------------
  ! Section 3 (bitmap definition section)
  ! integer numbers
  !--------------------------------------
  TYPE t_sec3
    SEQUENCE
    INTEGER :: bitmap  ! 0:internal bitmap 1: external bitmap number
    INTEGER :: missing ! missing value in integer field
  END TYPE t_sec3

  !--------------------------------------
  ! Section 3 (bitmap definition section)
  ! real numbers
  !--------------------------------------
  TYPE t_rsec3
    SEQUENCE
    REAL(dp) :: ignored
    REAL(dp) :: missing ! missing value in real field
  END TYPE t_rsec3

  !-----------------
  ! Section 4 (data)
  !-----------------
  TYPE t_sec4
    SEQUENCE
    INTEGER :: n_data      !  1 | Number of data values in array PSEC4.
    INTEGER :: bits        !  2 | Number of bits used for each encoded value.
    INTEGER :: grid_type   !  3 | 0: Grid point data, 128: Sph.harm.coefficients
    INTEGER :: packing     !  4 | 0: Simple packing,   64: Complex packing.
    INTEGER :: data_repr   !  5 | 0: Floating point,   32: Integer data.
    INTEGER :: flags       !  6 | 0: no                16: Additional flags.
    INTEGER :: reserved    !  7 | Reserved. Set to 0.
    INTEGER :: matrix      !  8 | 0: single datum      64: Matrix at each grid point.
    INTEGER :: bitmap2     !  9 | 0: no                32: Secondary bitmaps present.
    INTEGER :: width       ! 10 | 0: constant width    16: different width of second order values.
    INTEGER :: bits2       ! 11 | Number of bits for second order values.
    INTEGER :: wmo      (12: 15)! Reserved for WMO reserved flag fields. Set to 0.
    INTEGER :: start_cplx  ! 16 | For complex packing, start of packed data values.
    INTEGER :: scale_cplx  ! 17 | For complex packing, the scaling factor P.
    INTEGER :: j_cplx      ! 18 | For complex packing, pentagonal resolution parameter J.
    INTEGER :: k_cplx      ! 19 | For complex packing, pentagonal resolution parameter K.
    INTEGER :: m_cplx      ! 20 | For complex packing, pentagonal resolution parameter M.
    INTEGER :: non_miss    ! 21 | The number of non-missing values.
    INTEGER :: reserved2(22: 33)! Reserved. Set to 0.
    INTEGER :: o_coded     ! 34 | offset bit pointer to coded values in GRIB record. (returned by 'G', 'I' or 'J')
    INTEGER :: remaining(35:512)!
  END TYPE t_sec4

  !==============
  ! GRIB 1 record
  !==============
  INTEGER, PARAMETER :: INVALID_HANDLE = -1

  TYPE t_grib1
    !-----------------------------------------------
    ! parameters concerning the location in the file
    !-----------------------------------------------
    INTEGER            :: file         ! file handle
    INTEGER            :: krec         ! record number
#ifdef USE_PBSEEK64
    INTEGER(i8)        :: kpos         ! position in file
#else
    INTEGER            :: kpos         ! position in file
#endif
    INTEGER            :: blen         ! size of GRIB block in bytes
    !--------------
    ! Encoded data:
    !--------------
    INTEGER  ,POINTER  :: kgrib (:) => NULL()
    INTEGER            :: kleng     = 0 ! size (words) required for KGRIB
    INTEGER            :: kword     = 0 ! elements of KGRIB occupied by coded data.
    !----------------------------------------------
    ! Raw access to octets (get_octets, put_octets)
    !----------------------------------------------
    INTEGER ,POINTER   :: octets(:) => NULL()
    INTEGER ,POINTER   :: os0   (:) => NULL()
    INTEGER ,POINTER   :: os1   (:) => NULL()
    INTEGER ,POINTER   :: os2   (:) => NULL()
    INTEGER ,POINTER   :: os3   (:) => NULL()
    INTEGER ,POINTER   :: os4   (:) => NULL()
    INTEGER ,POINTER   :: os5   (:) => NULL()
    INTEGER ,POINTER   :: os6   (:) => NULL()
    !-------------------------------------------------------
    ! GRIB API handle (extension for migration to GRIB2 API)
    !-------------------------------------------------------
    INTEGER            :: handle    = INVALID_HANDLE
    LOGICAL            :: pbio      = .true.            ! Use PBIO
    !--------------
    ! Decoded data:
    !--------------
    CHARACTER          :: hoper     =' '! parameter used for decoding
    !----------
    ! section 0
    !----------
    TYPE (t_sec0)      :: isec0 ! (2)
    !----------------
    ! section 1 (PDB)
    !----------------
    TYPE (t_sec1)      :: isec1 ! (1024)
    TYPE (t_s1_dwd)    :: s1_dwd
    TYPE (t_s1_ecmwf)  :: s1_ecmwf
    INTEGER            :: pdtn          ! GRIB2 product definition template #
    INTEGER            :: dis           ! GRIB2 parameter discipline
    INTEGER            :: cat           ! GRIB2 parameter category
    INTEGER            :: num           ! GRIB2 parameter number
    INTEGER            :: ctyp          ! GRIB2 parameter constituentType
    !----------------
    ! section 2 (GDB)
    !----------------
    INTEGER            :: isec2  (4096) ! Integer data
    TYPE (t_rsec2)     :: rsec2 ! (512) ! Real data
    TYPE (t_s2_latlon) :: latlon        ! Integer data, lat/lon    grid
    TYPE (t_s2_gauss)  :: gauss         ! Integer data, Gauss      grid
    TYPE (t_s2_tri)    :: tri           ! Integer data, triangular grid
    TYPE (t_s2_sph)    :: sph           ! Integer data, sph. harm. grid
    !-------------------
    ! section 3 (bitmap)
    !-------------------
    TYPE (t_sec3)      :: isec3 ! (2)   ! bitmap
    TYPE (t_rsec3)     :: rsec3 ! (2)
    !-----------------
    ! section 4 (data)
    !-----------------
    TYPE (t_sec4)      :: isec4 !(512)        ! binary data
    REAL(dp) _POINTER  :: rsec4 (:) => NULL() ! decoded Real data
    INTEGER            :: klenp    = 0        ! The number of elements in array RSEC4
  END TYPE t_grib1

  !===========
  ! Interfaces
  !===========
  !------------------------------
  ! f90 wrapper for EMOS routines
  !------------------------------
  INTERFACE pbopen
    MODULE PROCEDURE pbopen_90
  END INTERFACE pbopen

  INTERFACE pbsetraw
    MODULE PROCEDURE pbsetraw_90
  END INTERFACE pbsetraw

  INTERFACE pbclose
    MODULE PROCEDURE pbclose_90
  END INTERFACE pbclose

  INTERFACE pbseek
#ifdef USE_PBSEEK64
    MODULE PROCEDURE pbseek64_90
#else
    MODULE PROCEDURE pbseek_90
#endif
  END INTERFACE pbseek

#ifdef USE_PBSEEK64
  INTERFACE pbtell
    MODULE PROCEDURE pbtell64
  END INTERFACE
#endif

  INTERFACE pbgrib
    MODULE PROCEDURE pbgrib_90
  END INTERFACE pbgrib

  INTERFACE pbwrite
    MODULE PROCEDURE pbwrite_90
  END INTERFACE pbwrite

  INTERFACE gribex
    MODULE PROCEDURE gribex_90
    MODULE PROCEDURE gribex_90_1
    MODULE PROCEDURE gribex_90_2
    MODULE PROCEDURE gribex_90_3
    MODULE PROCEDURE gribex_90_4
  END INTERFACE gribex

  !----------------------------------------
  ! Interfaces for public module procedures
  !----------------------------------------
  INTERFACE local_dwd
    MODULE PROCEDURE local_dwd_grib1
  END INTERFACE local_dwd

  INTERFACE local_ecmwf
    MODULE PROCEDURE local_ecmwf_grib1
  END INTERFACE local_ecmwf

  !-----------------------------------------
  ! Interfaces for private module procedures
  !-----------------------------------------
  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE   sec2_gauss
    MODULE PROCEDURE   sec2_latlon
    MODULE PROCEDURE   sec2_tri
    MODULE PROCEDURE   sec2_sph
    MODULE PROCEDURE  gauss_sec2
    MODULE PROCEDURE latlon_sec2
    MODULE PROCEDURE    tri_sec2
    MODULE PROCEDURE    sph_sec2
    MODULE PROCEDURE assign_grib1
  END INTERFACE !ASSIGNMENT (=)

  INTERFACE p_bcast
    MODULE PROCEDURE bcast_grib1
  END INTERFACE p_bcast

  INTERFACE destruct
    MODULE PROCEDURE destruct_grib
  END INTERFACE destruct

  !=========================
  ! private module variables
  !=========================
  integer :: size_grib = 0

  logical, parameter :: grib_verbose = .false.
! logical, parameter :: grib_verbose = .true.   ! Debug pbopen, gribex

!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE pbopen_90 (grib, filename, mode, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  CHARACTER(len=*)  ,INTENT(in)    :: filename
  CHARACTER(len=*)  ,INTENT(in)    :: mode     ! 'r','r+','w','a'
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL pbopen
    INTEGER :: iret
    if (grib_verbose) then
       write(0,*) "pbopen_90: file='",trim(filename),"', mode='",trim(mode),"'"
    end if
    CALL pbopen (grib% file, filename, mode, iret)
    grib% krec = 0
    grib% kpos = 0
    grib% pbio = .true.
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0)
        !                          the open was successful
      CASE (-1)
        CALL finish ('pbopen_90', 'could not open file "'//trim(filename)//&
                                  '", mode="'//trim (mode)//'"')
      CASE (-2)
        CALL finish ('pbopen_90', 'invalid filename: '//trim (filename))
      CASE (-3)
        CALL finish ('pbopen_90', 'invalid open mode specified: '//trim (mode))
      CASE default
        CALL finish ('pbopen_90', 'unexpected error condition')
      END SELECT
    ENDIF
#else
    call finish ("pbopen_90","cgribex not linked!")
#endif
  END SUBROUTINE pbopen_90
!------------------------------------------------------------------------------
  SUBROUTINE pbsetraw_90 (grib, raw)
  TYPE (t_grib1) ,INTENT(in) :: grib
  INTEGER        ,INTENT(in) :: raw
#ifndef NO_GRIBEX
    EXTERNAL pbsetraw
    if (.not. grib% pbio) call finish ("pbsetraw_90","file must be pbopen'ed")
    CALL pbsetraw (grib% file, raw)
#else
    call finish ("pbsetraw_90","cgribex not linked!")
#endif
  END SUBROUTINE pbsetraw_90
!------------------------------------------------------------------------------
  SUBROUTINE pbclose_90 (grib, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL pbclose
    INTEGER :: iret
    if (.not. grib% pbio) call finish ("pbclose_90","file must be pbopen'ed")
    CALL pbclose  (grib% file, iret)
    CALL destruct (grib)
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0)
        !                           the close was successful
      CASE (-1)
        CALL finish ('pbclose_90', 'an error occurred when closing the file')
      CASE default
        CALL finish ('pbclose_90', 'unexpected error condition')
      END SELECT
    ENDIF
#else
    call finish ("pbclose_90","cgribex not linked!")
#endif
  END SUBROUTINE pbclose_90
!------------------------------------------------------------------------------
  SUBROUTINE destruct_grib (grib)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
    grib% file = -1
    grib% krec = -1
    grib% kpos = -1
    IF (ASSOCIATED (grib% kgrib )) DEALLOCATE (grib% kgrib)
    IF (ASSOCIATED (grib% rsec4 )) DEALLOCATE (grib% rsec4)
    IF (ASSOCIATED (grib% octets)) DEALLOCATE (grib% octets)
    NULLIFY (grib% os0)
    NULLIFY (grib% os1)
    NULLIFY (grib% os2)
    NULLIFY (grib% os3)
    NULLIFY (grib% os4)
    NULLIFY (grib% os5)
  END SUBROUTINE destruct_grib
!------------------------------------------------------------------------------
#ifndef USE_PBSEEK64
  SUBROUTINE pbseek_90 (grib, koffset, kstart, kret)
  TYPE (t_grib1)    ,INTENT(in)    :: grib
  INTEGER           ,INTENT(in)    :: koffset
  INTEGER ,OPTIONAL ,INTENT(in)    :: kstart
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL pbseek
    INTEGER :: iret
    INTEGER :: istart
    if (.not. grib% pbio) call finish ("pbseek_90","file must be pbopen'ed")
    istart = 0; IF (PRESENT(kstart)) istart = kstart
    CALL pbseek (grib% file, koffset, istart, iret)
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0:)
        !                          byte offset from the start of the file
      CASE (-1)
        CALL finish ('pbseek_90', 'end-of-file is hit')
      CASE (-2)
        CALL finish ('pbseek_90', 'error in positioning the file')
      CASE default
        CALL finish ('pbseek_90', 'unexpected error condition')
      END SELECT
    ENDIF
#else
    call finish ("pbseek_90","cgribex not linked!")
#endif
  END SUBROUTINE pbseek_90
#endif
!------------------------------------------------------------------------------
#ifdef USE_PBSEEK64
  SUBROUTINE pbseek64_90 (grib, koffset, kstart, kret)
  TYPE (t_grib1)        ,INTENT(in)    :: grib
  INTEGER(i8)           ,INTENT(in)    :: koffset
  INTEGER     ,OPTIONAL ,INTENT(in)    :: kstart    ! == "whence" (0, 1 or 2)
  INTEGER(i8) ,OPTIONAL ,INTENT(out)   :: kret

    INTEGER(i8) :: iret
    INTEGER     :: whence
    EXTERNAL    :: pbseek64
    if (.not. grib% pbio) call finish ("pbseek64_90","file must be pbopen'ed")
    whence = 0; IF (PRESENT(kstart)) whence = kstart
    CALL pbseek64 (grib% file, koffset, whence, iret)
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0_i8:)
        !                          byte offset from the start of the file
      CASE (-1_i8)
        CALL finish ('pbseek64_90', 'end-of-file is hit')
      CASE (-2_i8)
        CALL finish ('pbseek64_90', 'error in positioning the file')
      CASE default
        CALL finish ('pbseek64_90', 'unexpected error condition')
      END SELECT
    ENDIF
  END SUBROUTINE pbseek64_90
#endif
!------------------------------------------------------------------------------
#ifdef USE_PBSEEK64
  subroutine pbtell64 (file, kpos)
    !------------------------------------------------
    ! Emulate 64-bit version of pbtell using pbseek64
    !------------------------------------------------
    integer,     intent(in)  :: file
    integer(i8), intent(out) :: kpos

    EXTERNAL           :: pbseek64
    integer, parameter :: SEEK_CUR = 1
    call pbseek64 (file, 0_i8, SEEK_CUR, kpos)
  end subroutine pbtell64
#endif
!------------------------------------------------------------------------------
  SUBROUTINE pbgrib_90 (grib, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL pbgrib
    INTEGER :: iret
    INTEGER :: koutlen
    !----------------------------
    ! determine size and position
    !----------------------------
    if (.not. grib% pbio) call finish ("pbgrib_90","file must be pbopen'ed")
    CALL PBSIZE (grib% file, grib% blen)
    if (grib% blen < 0) call finish ("pbgrib_90", "PBSIZE failed")
    CALL PBTELL (grib% file, grib% kpos)
    if (grib% kpos < 0) call finish ("pbgrib_90", "PBTELL failed")
    !--------------------------------------------
    ! (re)allocate buffer KGRIB with correct size
    !--------------------------------------------
    call reallocate_buffer (grib, (grib% blen+3)/4 )
    !-----------------
    ! read GRIB record
    !-----------------
    CALL PBGRIB (grib% file, grib% kgrib, grib% kleng*4, koutlen, iret )

    !+++++++++++++++++++++++++++++++++++++++++++
    ! work around wrong value returned by PBSIZE
    !+++++++++++++++++++++++++++++++++++++++++++
    if (iret==-3) then
      write (6,*) 'pbgrib_90: the size of KARRAY is not sufficient'
      write (6,*) 'kpos, blen, koutlen, iret =',grib% kpos, grib% blen, koutlen, iret
      write (6,*) 'RETRY:'
      CALL PBSEEK (grib, grib% kpos)
      grib% blen=koutlen
      call reallocate_buffer (grib, (grib% blen+3)/4 )
      CALL PBGRIB (grib% file, grib% kgrib, grib% kleng*4, koutlen, iret )
      write (6,*) 'kpos, blen, koutlen, iret =',grib% kpos, grib% blen, koutlen, iret
    endif

    grib% kword = (koutlen + 3) / 4
    IF (iret == 0) grib% krec = grib% krec + 1          ! Bump record number
    !-----------------
    ! error processing
    !-----------------
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0)
        !                          GRIB product has been read successfully
      CASE (-1)
        CALL finish ('pbgrib_90', 'end-of-file was hit')
      CASE (-2)
        CALL finish ('pbgrib_90', 'error in the file-handling')
      CASE (-3)
        CALL finish ('pbgrib_90', 'the size of KARRAY is not sufficient')
      CASE default
        CALL finish ('pbgrib_90', 'unexpected error condition')
      END SELECT
    ENDIF
#else
    call finish ("pbgrib_90","cgribex not linked!")
#endif
  END SUBROUTINE pbgrib_90
!------------------------------------------------------------------------------
  subroutine pbwrite_90 (grib, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL pbwrite
    INTEGER :: iret
    INTEGER :: koutlen
    if (.not. grib% pbio) call finish ("pbwrite_90","file must be pbopen'ed")
    koutlen = grib% blen
    call pbwrite (grib% file, grib% kgrib, koutlen, iret)
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      SELECT CASE (iret)
      CASE (0:)
        !                          number of bytes written to the file
      CASE (-1)
        CALL finish ('pbwrite_90', 'error writing to the file')
      CASE default
        CALL finish ('pbwrite_90', 'unexpected error condition')
      END SELECT
    ENDIF
#else
    call finish ("pbwrite_90","cgribex not linked!")
#endif
  end subroutine pbwrite_90
!------------------------------------------------------------------------------
  SUBROUTINE gribex_90_1 (grib, hoper, x, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  CHARACTER         ,INTENT(in)    :: hoper
  REAL(wp)          ,INTENT(inout) :: x (:)
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
    INTEGER :: n
    if (.not. grib% pbio) call finish ("gribex_90_1","file must be pbopen'ed")
    SELECT CASE (hoper)
    CASE ('C')
      n = SIZE(x)
      CALL reallocate_data (grib,n)
      grib% rsec4 (1:n) = RESHAPE (x, (/n/))
      CALL gribex_90 (grib, hoper, kret)
    END SELECT
    SELECT CASE (hoper)
    CASE ('D','R')
      CALL gribex_90 (grib, hoper, kret)
      n = grib% isec4% n_data
      IF (n /= SIZE(x)) THEN
        WRITE (0,*) 'GRIB record, position:', grib% krec, grib% kpos
        WRITE (0,*) 'isec4(1) =',grib% isec4% n_data, &
          '/= size(x) =',SIZE(x),', shape(x) =',SHAPE(x)
        CALL finish ('gribex_90_1 ('//hoper//')', 'data size mismatch?')
      ENDIF
      x = RESHAPE (grib% rsec4 (1:n), SHAPE (x))
    END SELECT
  END SUBROUTINE gribex_90_1
!------------------------------------------------------------------------------
  SUBROUTINE gribex_90_2 (grib, hoper, x, kret, scanmode)
  TYPE (t_grib1)    ,INTENT(inout) :: grib      ! GRIB derived data type
  CHARACTER         ,INTENT(in)    :: hoper     ! 'C'ode, 'D'ecode, 'R'
  REAL(wp)          ,INTENT(inout) :: x (:,:)   ! field to extract from GRIB
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret      ! error return value
  INTEGER ,OPTIONAL ,INTENT(in)    :: scanmode  ! scanning mode for x
    INTEGER :: n, sm
    if (.not. grib% pbio) call finish ("gribex_90_2","file must be pbopen'ed")
    SELECT CASE (hoper)
    CASE ('C')
      if(present(scanmode)) call finish ('gribex_90_2','scanmode + hoper=C')
      n = SIZE(x)
      CALL reallocate_data (grib,n)
      grib% rsec4 (1:n) = RESHAPE (x, (/n/))
      CALL gribex_90 (grib, hoper, kret)
    END SELECT
    SELECT CASE (hoper)
    CASE ('D','R')
      CALL gribex_90 (grib, hoper, kret)
      n = grib% isec4% n_data
      IF (n /= SIZE(x)) THEN
        WRITE (0,*) 'GRIB record, position:', grib% krec, grib% kpos
        WRITE (0,*) 'isec4(1) =',grib% isec4% n_data, &
          '/= size(x) =',SIZE(x),', shape(x) =',SHAPE(x)
        CALL finish ('gribex_90_2 ('//hoper//')', 'data size mismatch?')
      ENDIF
      x = RESHAPE (grib% rsec4 (1:n), SHAPE (x))
      if (present(scanmode)) then
        select case (grib% isec2(1))
        case (WMO6_LATLON, WMO6_GAUSSIAN, WMO6_ROTLL)
          sm  = ieor(scanmode,grib% isec2(11))
          if (iand(scanmode, WMO8_J_CONSECUTIVE) /= 0) &
            call finish ('gribex_90_2','scanmode==WMO8_J_CONSECUTIVE')
          if (iand(sm,       WMO8_J_CONSECUTIVE) /= 0) &
            call finish ('gribex_90_2','scan_mode==WMO8_J_CONSECUTIVE')
          if (iand(sm,       WMO8_I_NEGATIVE) /= 0) x = x(  size(x,1):1:-1,:)
          if (iand(sm,       WMO8_J_POSITIVE) /= 0) x = x(:,size(x,2):1:-1  )
        case default
          if (grib_verbose) &
             call message ('gribex_90_2','scanmode without LATLON or GAUSSIAN')
!         call finish ('gribex_90_2','scanmode without LATLON or GAUSSIAN')
        end select
      endif
    END SELECT
  END SUBROUTINE gribex_90_2
!------------------------------------------------------------------------------
  SUBROUTINE gribex_90_3 (grib, hoper, x, kret, scanmode)
  TYPE (t_grib1)    ,INTENT(inout) :: grib      ! GRIB derived data type
  CHARACTER         ,INTENT(in)    :: hoper     ! 'C'ode, 'D'ecode, 'R'
  REAL(wp)          ,INTENT(inout) :: x (:,:,:) ! field to extract from GRIB
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret      ! error return value
  INTEGER ,OPTIONAL ,INTENT(in)    :: scanmode  ! scanning mode for x
    INTEGER :: n, sm
    if (.not. grib% pbio) call finish ("gribex_90_3","file must be pbopen'ed")
    SELECT CASE (hoper)
    CASE ('C')
      if(present(scanmode)) call finish ('gribex_90_3','scanmode + hoper=C')
      n = SIZE(x)
      CALL reallocate_data (grib,n)
      grib% rsec4 (1:n) = RESHAPE (x, (/n/))
      CALL gribex_90 (grib, hoper, kret)
    END SELECT
    SELECT CASE (hoper)
    CASE ('D','R')
      CALL gribex_90 (grib, hoper, kret)
      n = grib% isec4% n_data
      IF (n /= SIZE(x)) THEN
        WRITE (0,*) 'GRIB record, position:', grib% krec, grib% kpos
        WRITE (0,*) 'isec4(1) =',grib% isec4% n_data, &
          '/= size(x) =',SIZE(x),', shape(x) =',SHAPE(x)
        WRITE (6,*) 'GRIB record, position:', grib% krec, grib% kpos
        WRITE (6,*) 'klenp    =',grib% klenp
        WRITE (6,*) 'isec4(1) =',grib% isec4% n_data, &
          '/= size(x) =',SIZE(x),', shape(x) =',SHAPE(x)
        call print_isec4 (grib% isec4)
        CALL finish ('gribex_90_3 ('//hoper//')', 'data size mismatch?')
      ENDIF
      x = RESHAPE (grib% rsec4 (1:n), SHAPE (x))
      if (present(scanmode)) then
        select case (grib% isec2(1))
        case (WMO6_LATLON, WMO6_GAUSSIAN, WMO6_ROTLL)
          sm  = ieor(scanmode,grib% isec2(11))
          if (iand(scanmode, WMO8_J_CONSECUTIVE) /= 0) &
            call finish ('gribex_90_3','scanmode==WMO8_J_CONSECUTIVE')
          if (iand(sm,       WMO8_J_CONSECUTIVE) /= 0) &
            call finish ('gribex_90_3','scan_mode==WMO8_J_CONSECUTIVE')
          if (iand(sm,       WMO8_I_NEGATIVE) /= 0) x = x(  size(x,1):1:-1,:,:)
          if (iand(sm,       WMO8_J_POSITIVE) /= 0) x = x(:,size(x,2):1:-1,:  )
        case default
          if (grib_verbose) &
             call message ('gribex_90_3','scanmode without LATLON or GAUSSIAN')
!         call finish ('gribex_90_3','scanmode without LATLON or GAUSSIAN')
        end select
      endif
    END SELECT
  END SUBROUTINE gribex_90_3
!------------------------------------------------------------------------------
  SUBROUTINE gribex_90_4 (grib, hoper, x, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  CHARACTER         ,INTENT(in)    :: hoper
  REAL(wp)          ,INTENT(inout) :: x (:,:,:,:)
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
    INTEGER :: n
    if (.not. grib% pbio) call finish ("gribex_90_4","file must be pbopen'ed")
    SELECT CASE (hoper)
    CASE ('C')
      n = SIZE(x)
      CALL reallocate_data (grib,n)
      grib% rsec4 (1:n) = RESHAPE (x, (/n/))
      CALL gribex_90 (grib, hoper, kret)
    END SELECT
    SELECT CASE (hoper)
    CASE ('D','R')
      CALL gribex_90 (grib, hoper, kret)
      n = grib% isec4% n_data
      IF (n /= SIZE(x)) THEN
        WRITE (0,*) 'GRIB record, position:', grib% krec, grib% kpos
        WRITE (0,*) 'isec4(1) =',grib% isec4% n_data, &
          '/= size(x) =',SIZE(x),', shape(x) =',SHAPE(x)
        CALL finish ('gribex_90_4 ('//hoper//')', 'data size mismatch?')
      ENDIF
      x = RESHAPE (grib% rsec4 (1:n), SHAPE (x))
    END SELECT
  END SUBROUTINE gribex_90_4
!------------------------------------------------------------------------------
  SUBROUTINE gribex_90 (grib, hoper, kret)
  TYPE (t_grib1)    ,INTENT(inout) :: grib
  CHARACTER         ,INTENT(in)    :: hoper
  INTEGER ,OPTIONAL ,INTENT(out)   :: kret
#ifndef NO_GRIBEX
    EXTERNAL gribex
    target           :: grib
    INTEGER          :: iret, ntry, nx, n_data, reduced
    integer          :: l_enc, l_dec
    integer ,pointer :: lu(:)
    integer          :: nsec(0:5)       ! Debug: sizes of sections
    external         :: grprs0, grprs1, grprs2
!   integer          :: version
    integer          :: s, a, b  ! decode reference value
    real(sp)         :: r        ! GRIB reference value
    logical          :: allsame  ! flag for all data of same value

    if (.not. grib% pbio) call finish ("gribex_90","file must be pbopen'ed")
    if (hoper == 'V') then
       !----------------------------------------------------------------------
       ! Print GRIBEX version and return (needs sufficiently new grib library)
       !----------------------------------------------------------------------
       CALL GRIBEX (grib% isec0, &! PINT
                    grib% isec1, &! PINT
                    grib% isec2, &! PINT
                    grib% rsec2, &! PVOID
                    grib% isec3, &! PINT
                    grib% rsec3, &! PVOID
                    grib% isec4, &! PINT
                    grib% rsec4, &! PVOID
                    grib% klenp, &! INT
                    grib% kgrib, &! PINT
                    grib% kleng, &! INT
                    grib% kword, &! PINT
                    'V'        , &! STRING
                    iret)         ! PINT
       if (present (kret)) kret = iret
       return
    end if
    !------------------------------------------
    ! copy integer data to grid specific buffer
    !------------------------------------------
    SELECT CASE (hoper)
    CASE ('C')
      SELECT CASE (grib% isec2 (1))
      CASE (WMO6_LATLON, WMO6_ROTLL)
        grib% isec2 = grib% latlon
      CASE (WMO6_GAUSSIAN)
        grib% isec2 = grib% gauss
      CASE (WMO6_HARMONIC)
        grib% isec2 = grib% sph
      CASE (DWD6_ICOSAHEDRON)
        grib% isec2 = grib% tri
      END SELECT
      if (local_dwd(grib).and.                                 &
          grib% s1_dwd% local_ident == grib% isec1% local_ident) then
        lu => grib% isec1% local_use
        grib% isec1% local_ident = grib% s1_dwd% local_ident
        lu(38)  =       grib% s1_dwd% day_number
        lu(39)  =       grib% s1_dwd% record_number / 65536
        lu(40)  =  mod (grib% s1_dwd% record_number /   256, 256)
        lu(41)  =  mod (grib% s1_dwd% record_number        , 256)
        lu(42)  =       grib% s1_dwd% decoding
        lu(43)  =       grib% s1_dwd% element_no
        lu(44)  =       grib% s1_dwd% year
        lu(45)  =       grib% s1_dwd% month
        lu(46)  =       grib% s1_dwd% day
        lu(47)  =       grib% s1_dwd% hour
        lu(48)  =       grib% s1_dwd% minute
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
! fix bug in MPI GRIB library
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
!       lu(49)  =       grib% s1_dwd% exp
!       lu(50)  =       grib% s1_dwd% run_type
        lu(49)  =       grib% s1_dwd% exp            &
                +       grib% s1_dwd% run_type * 2**14
        lu(50)  =       0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
        lu(51)  =       grib% s1_dwd% user_id
        lu(52)  =       grib% s1_dwd% experiment_id
        lu(53)  =       grib% s1_dwd% ensemble_id
        lu(54)  =       grib% s1_dwd% ensemble_size
        lu(55)  =       grib% s1_dwd% ensemble_no
        lu(56)  =       grib% s1_dwd% major_version
        lu(57)  =       grib% s1_dwd% minor_version
        lu(58)  =       0           ! blank (even block length)
      endif

      if (local_ecmwf (grib).and.                                 &
          grib% s1_ecmwf% local_ident == grib% isec1% local_ident) then
        grib% isec1% reserved(:) = 0
        call encode_ecmwf (grib% s1_ecmwf, &
                           grib% isec1% local_use)
      endif
    END SELECT
    reduced = 0
    ntry    = 0
    allsame = .false.
    DO
      !--------------------------------------
      ! (re)allocate data buffer for decoding
      !--------------------------------------
      SELECT CASE (hoper)
      CASE ('D','R')
        !++++++++++++++++++++++
        ! decode metadata first
        !++++++++++++++++++++++
        iret = 0
        CALL reallocate_data (grib, 0)
FTRACE_BEGIN("gribex_90:GRIBEX:J")
        CALL GRIBEX (grib% isec0, &! PINT
                     grib% isec1, &! PINT
                     grib% isec2, &! PINT
                     grib% rsec2, &! PVOID
                     grib% isec3, &! PINT
                     grib% rsec3, &! PVOID
                     grib% isec4, &! PINT
                     grib% rsec4, &! PVOID
                     grib% klenp, &! INT
                     grib% kgrib, &! PINT
                     grib% kleng, &! INT
                     grib% kword, &! PINT
                     'J'        , &! STRING
                     iret)         ! PINT
FTRACE_END  ("gribex_90:GRIBEX:J")
        if (iret /= 0) call finish ("gribex_90","decoding of metadata failed")
        n_data  = grib% isec4% n_data
        reduced = grib% isec2(17)
        nx      = grib% isec2(2)
        if (nx < 0 .or. nx > 2**30) then
           write(0,*) "gribex_90: nx      =", nx
           write(0,*) "gribex_90: reduced =", reduced
           write(0,*) "gribex_90: n_data  =", n_data
           call grprs0 (grib% isec0)
           call grprs1 (grib% isec0, grib% isec1)
           call grprs2 (grib% isec0, grib% isec2, grib% rsec2)
           print *
           print *, "Section 2, elements 1..22:"
           print *, grib% isec2(1:22)
           print *
           l_dec = bufr_size_enc (grib, verbose=2, nsec=nsec)
           print *
           write(6,'(a,i6)') ' Size of GRIB message:',l_dec
           write(6,'(a,i6)') ' Size of section 0   :',nsec(0)
           write(6,'(a,i6)') ' Size of section 1   :',nsec(1)
           write(6,'(a,i6)') ' Size of section 2   :',nsec(2)
           write(6,'(a,i6)') ' Size of section 3   :',nsec(3)
           write(6,'(a,i6)') ' Size of section 4   :',nsec(4)
           write(6,'(a,i6)') ' Size of section 5   :',nsec(5)
           call finish ("gribex_90","bad grib or bug in GRIB library?")
        end if
        if (reduced == 1 .and. hoper == 'R') then
           !-------------------------------------
           ! Reduced grid: prepare to expand grid
           !-------------------------------------
           nx = grib% isec2(2)
           if (nx == 0 .or. grib% isec2(6) == 0) then
              nx = maxval (grib% isec2(23:22+grib% isec2(3)))
              call message('gribex_90 ('//hoper//')', &
                   ' reduced grid, using maxval (pl) as number of longitudes')
              if (grib_verbose) then
                 write(0,*) "gribex_90: ny =", grib% isec2(3)
                 write(0,*) "gribex_90: nx =", nx
              end if
           end if
           n_data = nx * grib% isec2(3)
        end if
        if (reduced == 1 .and. hoper /= 'R') then
           call message('gribex_90','reduced grid, but HOPER = '//hoper)
           if (grib% isec2(2) * grib% isec2(3) /= n_data) then
              write (0,*) 'gribex_90: section 2: ni    =', grib% isec2(2)
              write (0,*) 'gribex_90: section 2: nj    =', grib% isec2(3)
              write (0,*) 'gribex_90: section 2: resol =', grib% isec2(6)
              write (0,*) 'gribex_90: isec4% n_data    =', n_data
           end if
        end if
        if (n_data == 0 .and. grib% isec4% bits     == 0 &
                        .and. grib% isec4% non_miss == 0) then
          !-------------------------------------------------------------
          ! all values have the same value, deduce n_data from section 2
          !-------------------------------------------------------------
          allsame = .true.
          select case (grib% isec2 (1))
          case (WMO6_LATLON, WMO6_ROTLL)
            n_data =  grib% latlon% ni    * grib% latlon% nj
          case (WMO6_GAUSSIAN)
            n_data =  grib% gauss%  ni    * grib% gauss%  nj
          case (DWD6_ICOSAHEDRON)
            n_data = (grib% tri% ni+1)**2 * grib% tri%    nd
          end select
          write(6,*) 'all data has the same value, n_data =',n_data
        endif
        CALL reallocate_data (grib, n_data)
        grib% hoper            = hoper
        grib% isec1% local_use = -huge(0)
      CASE ('C')
        l_dec = bufr_size_dec (grib,0)
        if (grib_verbose) then
           l_dec = bufr_size_dec (grib, verbose=2)
        end if
        grib% blen = l_dec
!       CALL reallocate_buffer (grib, (l_dec+3)/4)
        CALL reallocate_buffer (grib, l_dec/4+1  +1)!???
      CASE default
        CALL reallocate_data   (grib, 0)
        CALL reallocate_buffer (grib, 0)
      END SELECT
      !----------------
      ! decode / encode
      !----------------
      iret = 1
!     !++++++++++++++++++++++++++++++++++++++++++++
!     ! error handling of libgrib does not work yet
!     !++++++++++++++++++++++++++++++++++++++++++++
!     iret = 0
FTRACE_BEGIN("gribex_90:GRIBEX:"//hoper)
      CALL GRIBEX (grib% isec0, &! PINT
                   grib% isec1, &! PINT
                   grib% isec2, &! PINT
                   grib% rsec2, &! PVOID
                   grib% isec3, &! PINT
                   grib% rsec3, &! PVOID
                   grib% isec4, &! PINT
                   grib% rsec4, &! PVOID
                   grib% klenp, &! INT
                   grib% kgrib, &! PINT
                   grib% kleng, &! INT
                   grib% kword, &! PINT
                   hoper,       &! STRING
                   iret)         ! PINT
FTRACE_END  ("gribex_90:GRIBEX:"//hoper)

      if (allsame) then
        !--------------------------------------------------
        ! all values have the same value (=reference value)
        !--------------------------------------------------
        call get_octets (grib)
        s = 1; if (iand (128, grib% os4(7)) /= 0) s = -1
        a =        iand (127, grib% os4(7))
        b = grib% os4(8)*256**2+grib% os4(9)*256+grib% os4(10)
        r = s * 2._sp**(-24) * b * 16._sp**(a-64)
!       write(6,*) 'all values have the same (reference) value'
!       write(6,*) 's =',s   ! sign
!       write(6,*) 'a =',a   ! exponent
!       write(6,*) 'b =',b   ! mantissa
!       write(6,*) 'r =',r   ! reference value
!       write(6,*) 'd =',grib% isec1% factor
        write(6,*) 'all data set to the reference value =',r
        grib% rsec4 = r
        grib% isec4% n_data = n_data
        deallocate (grib% octets)
      endif
      if (iret /= 0) then
        write (0,*)  'gribex_90('//hoper//'):  KRET from GRIBEX is ',iret
        call message('gribex_90('//hoper//')','KRET from GRIBEX is /= 0')
      endif
      !-------------------------------------------
      ! crosscheck sizes of encoded bufr
      ! and 3dvar assumption based on decoded bufr
      !-------------------------------------------
      l_dec = bufr_size_dec (grib, verbose=0)
      l_enc = bufr_size_enc (grib, verbose=0, nsec=nsec)
      if (l_dec/=l_enc) then
        if (l_enc == l_dec + 1                                     .and. &
            ((grib% isec1% local_ident == 253 .and. nsec(1) == 66) .or. &
             (grib% isec1% local_ident == 153 .and. nsec(1) == 70)      )) then
          !------------------------------------------------------------
          ! allow for difference of 1 octed in dwd pds1 local extension
          ! blank octet not written by MPIfM cgribex or DWD grib_api,
          ! but even number of octets required by WMO regulation 92.1.4
          !------------------------------------------------------------
        else
          l_dec = bufr_size_dec (grib, verbose=2)
          l_enc = bufr_size_enc (grib, verbose=2)
          write(0,*)
          write(0,*) "l_dec, l_enc =", l_dec, l_enc
          call finish ('gribex_90 ('//hoper//')','buffer size mismatch')
        endif
      endif
      !------------------------------
      ! retry if buffer was too small
      !------------------------------
      SELECT CASE (hoper)
      CASE ('D','R')
        IF (iret > 0 .AND. SIZE (grib% kgrib) < grib% isec4% n_data) then
           ntry = ntry + 1
           if (ntry >= 3) &
                call finish ('gribex_90 ('//hoper//')',"failing after 3 tries")
           CYCLE
        end IF
      END SELECT
      EXIT
    END DO
    !-----------------
    ! error processing
    !-----------------
    IF (PRESENT(kret)) THEN
      kret = iret
    ELSE
      IF (iret /= 0) CALL finish ('gribex_90',&
        'error when decoding grib record')
    ENDIF
    !------------------------------------------
    ! copy integer data to grid specific buffer
    !------------------------------------------
    SELECT CASE (hoper)
    CASE ('D','R','J')
      SELECT CASE (grib% isec2 (1))
      CASE (WMO6_LATLON, WMO6_ROTLL)
        grib%   latlon = grib% isec2
      CASE (WMO6_GAUSSIAN)
        grib%   gauss  = grib% isec2
      CASE (WMO6_HARMONIC)
        grib%   sph    = grib% isec2
      CASE (DWD6_ICOSAHEDRON)
        grib%   tri    = grib% isec2
      END SELECT
      if (local_dwd(grib)) then
        lu => grib% isec1% local_use
        grib% s1_dwd% local_ident   = grib% isec1% local_ident
        grib% s1_dwd% day_number    = lu(38)
        grib% s1_dwd% record_number = lu(39) * 65536 + &
                                      lu(40) *   256 + &
                                      lu(41)
        grib% s1_dwd% decoding      = lu(42)
        grib% s1_dwd% element_no    = lu(43)
        grib% s1_dwd% year          = lu(44)
        grib% s1_dwd% month         = lu(45)
        grib% s1_dwd% day           = lu(46)
        grib% s1_dwd% hour          = lu(47)
        grib% s1_dwd% minute        = lu(48)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
! fix bug in MPI GRIB library
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
!       grib% s1_dwd% exp           = lu(49)
!       grib% s1_dwd% run_type      = lu(50)
        grib% s1_dwd% run_type      = lu(50) * 2     &
                                    + lu(49) / 2**14
        grib% s1_dwd% exp     = mod ( lu(49) , 2**14 )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (grib% s1_dwd% local_ident == 253) then
        grib% s1_dwd% user_id       = lu(51)
        grib% s1_dwd% experiment_id = lu(52)
        grib% s1_dwd% ensemble_id   = lu(53)
        grib% s1_dwd% ensemble_size = lu(54)
        grib% s1_dwd% ensemble_no   = lu(55)
        grib% s1_dwd% major_version = lu(56)
        grib% s1_dwd% minor_version = lu(57)
        else
          grib% s1_dwd% user_id       = -1
          grib% s1_dwd% experiment_id = -1
          grib% s1_dwd% ensemble_id   = -1
          grib% s1_dwd% ensemble_size = -1
          grib% s1_dwd% ensemble_no   = -1
          grib% s1_dwd% major_version = -1
          grib% s1_dwd% minor_version = -1
        endif
      else
        grib% s1_dwd% local_ident   = -1
      endif
      if (local_ecmwf (grib)) then
        call decode_ecmwf (grib% s1_ecmwf, &
                           grib% isec1% local_ident, &
                           grib% isec1% local_use)
      else
        grib% s1_ecmwf% local_ident = -1
      endif
    END SELECT
#else
    call finish ("gribex_90","cgribex not linked!")
#endif
  END SUBROUTINE gribex_90
!==============================================================================
  subroutine grprs1lu (grib)
  TYPE (t_grib1) ,INTENT(in) ,target :: grib
  !========================================
  ! print local use part of section 1 (PDB)
  !========================================

!   integer ::i
    type (t_s1_dwd)  ,pointer :: dwd

    write(6,"()")
    write(6,"(a)")    ' Section 1 - Product Definition Section, Local Use.'
    write(6,"(a)")    ' --------------------------------------------------'
    write(6,"(a44,i3)") ' Originating centre. ', grib% isec1% center
    write(6,"(a44,i3)") ' Sub-centre          ', grib% isec1% sub_center
    write(6,"(a44,i3)") ' local_flag.         ', grib% isec1% local_flag
    write(6,"(a44,i3)") ' local_ident.        ', grib% isec1% local_ident

!  do i=lbound(grib% isec1% local_use,1), &
!       ubound(grib% isec1% local_use,1)
!    if (grib% isec1% local_use(i) /= -huge(0))             &
!      write(6,"(a,i2,a,i3)")                               &
!        ' local_use (',i,')                             ', &
!        grib% isec1% local_use(i)
!  end do

    dwd => grib% s1_dwd
    if (local_dwd(grib) .and.                                &
        grib% s1_dwd% local_ident == grib% isec1% local_ident) then
      write(6,"()")
      write(6,"(a39,i8)")'TAGESNUMMER  BACKUP-FILE     ', dwd% day_number
      write(6,"(a39,i8)")'RECORDNUMMER BACKUP-FILE     ', dwd% record_number
      write(6,"(a39,i8)")'ENTSCHLUESSELN MFA(1)/MFB(0) ', dwd% decoding
      write(6,"(a39,i8)")'ZUSATZ-ELEMENTNUMMER         ', dwd% element_no
      write(6,"(a39,i8)")'JAHR                         ', dwd% year
      write(6,"(a39,i8)")'MONAT                        ', dwd% month
      write(6,"(a39,i8)")'TAG                          ', dwd% day
      write(6,"(a39,i8)")'STUNDE                       ', dwd% hour
      write(6,"(a39,i8)")'MINUTE                       ', dwd% minute
      write(6,"(a39,i8)")'EXPERIMENT                   ', dwd% exp
      write(6,"(a39,i8)")'RUN (0=haupt, 2=ass, 3=test) ', dwd% run_type

      if (grib% s1_dwd% local_ident == 253) then
      write(6,"(a39,i8)")'Experiment identifier        ', dwd% experiment_id
      write(6,"(a39,i8)")'Ensemble identification      ', dwd% ensemble_id
      write(6,"(a39,i8)")'Number of ensemble members   ', dwd% ensemble_size
      write(6,"(a39,i8)")'Esemble member number        ', dwd% ensemble_no
      write(6,"(a39,i8)")'Model major version number   ', dwd% major_version
      write(6,"(a39,i8)")'Model minor version number   ', dwd% minor_version
      endif

    end if

    if (local_ecmwf(grib) .and. &
      grib% s1_ecmwf% local_ident == grib% isec1% local_ident) &
      call grprs1ec (grib% s1_ecmwf)

  end subroutine grprs1lu
!==============================================================================
  FUNCTION no_spectral_coefficients(nj, nk, nm) RESULT (nsp)
  INTEGER, INTENT(in) :: nj, nk, nm
  INTEGER             :: nsp

    INTEGER :: n, m
    INTEGER :: ns

    nsp = 0
    DO m = 0, nm
      IF (nj+m <= nk) THEN
        ns = nj+m
      ELSE
        ns = nk
      END IF
      DO n = m, ns
        nsp = nsp+1
      END DO
    END DO

  END FUNCTION no_spectral_coefficients
!==============================================================================
  SUBROUTINE sec2_gauss (y, x)
    INTEGER           ,INTENT(inout) :: y (:)
    TYPE (t_s2_gauss) ,INTENT(in)    :: x
      y (1:22) = TRANSFER (x, y(1:22))
  END SUBROUTINE sec2_gauss
!------------------------------------------------------------------------------
  SUBROUTINE gauss_sec2 (y, x)
    TYPE (t_s2_gauss) ,INTENT(out) :: y
    INTEGER           ,INTENT(in)  :: x (:)
      y = TRANSFER (x(1:22), y)
  END SUBROUTINE gauss_sec2
!------------------------------------------------------------------------------
  SUBROUTINE sec2_latlon (y, x)
    INTEGER            ,INTENT(inout) :: y (:)
    TYPE (t_s2_latlon) ,INTENT(in)    :: x
      y (1:22) = TRANSFER (x, y(1:22))
  END SUBROUTINE sec2_latlon
!------------------------------------------------------------------------------
  SUBROUTINE latlon_sec2 (y, x)
    TYPE (t_s2_latlon) ,INTENT(out) :: y
    INTEGER            ,INTENT(in)  :: x (:)
      y = TRANSFER (x(1:22), y)
  END SUBROUTINE latlon_sec2
!------------------------------------------------------------------------------
  SUBROUTINE sec2_tri (y, x)
    INTEGER         ,INTENT(inout) :: y (:)
    TYPE (t_s2_tri) ,INTENT(in)    :: x
      y (1:22) = TRANSFER (x, y(1:22))
  END SUBROUTINE sec2_tri
!------------------------------------------------------------------------------
  SUBROUTINE tri_sec2 (y, x)
    TYPE (t_s2_tri) ,INTENT(out) :: y
    INTEGER         ,INTENT(in)  :: x (:)
      y = TRANSFER (x(1:22), y)
  END SUBROUTINE tri_sec2
!------------------------------------------------------------------------------
  SUBROUTINE sec2_sph (y, x)
    INTEGER         ,INTENT(inout) :: y (:)
    TYPE (t_s2_sph) ,INTENT(in)    :: x
      y (1:22) = TRANSFER (x, y(1:22))
  END SUBROUTINE sec2_sph
!------------------------------------------------------------------------------
  SUBROUTINE sph_sec2 (y, x)
    TYPE (t_s2_sph) ,INTENT(out) :: y
    INTEGER         ,INTENT(in)  :: x (:)
      y = TRANSFER (x(1:22), y)
  END SUBROUTINE sph_sec2
!------------------------------------------------------------------------------
  !-------------------------------------------
  ! allocate 'kgrib' with the appropriate size
  ! preset   'kgrib' with zeros
  ! store size of 'kgrib' in 'kleng'
  !-------------------------------------------
  SUBROUTINE reallocate_buffer (grib, n)
  TYPE (t_grib1) ,INTENT(inout) :: grib
  INTEGER        ,INTENT(in)    :: n
    grib% kleng = n
    IF (.NOT. ASSOCIATED (grib% kgrib)) THEN
       ALLOCATE (grib% kgrib (grib% kleng))
       grib% kgrib = 0
    ELSE IF (SIZE (grib% kgrib) < grib% kleng) THEN
       DEALLOCATE (grib% kgrib)
       ALLOCATE   (grib% kgrib (grib% kleng))
       grib% kgrib = 0
    ENDIF
    grib% kleng = SIZE (grib% kgrib)
  END SUBROUTINE reallocate_buffer
!------------------------------------------------------------------------------
  !-------------------------------------------
  ! allocate 'rsec4' with the appropriate size
  ! preset   'rsec4' with zeros
  ! store size of 'rsec4' in 'klenp'
  !-------------------------------------------
  SUBROUTINE reallocate_data (grib, n)
  TYPE (t_grib1) ,INTENT(inout)  :: grib
  INTEGER        ,INTENT(in)     :: n
    grib% klenp = n
    IF (.NOT. ASSOCIATED (grib% rsec4)) then
       ALLOCATE (grib% rsec4 (grib% klenp))
       grib% rsec4 = 0._wp
    ELSE IF (SIZE (grib% rsec4) /= grib% klenp) THEN
       DEALLOCATE (grib% rsec4)
       ALLOCATE (grib% rsec4 (grib% klenp))
       grib% rsec4 = 0._wp
    ENDIF
    grib% klenp = SIZE (grib% rsec4)
  END SUBROUTINE reallocate_data
!==============================================================================
  !-------------------
  ! copy a GRIB record
  !-------------------
  SUBROUTINE assign_grib1 (y, x)
  TYPE (t_grib1) ,INTENT(inout)  :: y
  TYPE (t_grib1) ,INTENT(in)     :: x
    !-----------------------------------------------
    ! parameters concerning the location in the file
    !-----------------------------------------------
    !  file = unchanged             ! file index
    !  krec = unchanged             ! record number
    !  kpos = unchanged             ! position in file
    y% blen = x% blen               ! size of GRIB block in bytes
    !--------------
    ! Encoded data:
    !--------------
    y% kleng = x% kleng             ! size (words) required for kgrib
    y% kword = x% kword             ! elements of KGRIB occupied by data.
    if (associated (x% kgrib)) then ! encoded data
       call reallocate_buffer (y, size (x% kgrib))
       y% kgrib (:size (x% kgrib)) = x% kgrib(:)
    endif
    !--------------
    ! Decoded data:
    !--------------
    y% hoper = x% hoper
    !----------
    ! section 0
    !----------
    y% isec0 = x% isec0
    !----------------
    ! section 1 (PDB)
    !----------------
    y% isec1    = x% isec1
    y% s1_dwd   = x% s1_dwd
    y% s1_ecmwf = x% s1_ecmwf
    !----------------
    ! section 2 (GDB)
    !----------------
    y% isec2  = x% isec2  ! Integer data
    y% rsec2  = x% rsec2  ! Real data
    y% latlon = x% latlon ! Integer data, lat/lon    grid
    y% gauss  = x% gauss  ! Integer data, Gauss      grid
    y% tri    = x% tri    ! Integer data, triangular grid
    y% sph    = x% sph    ! Integer data, sph. harm. grid
    !-------------------
    ! section 3 (bitmap)
    !-------------------
    y% isec3 = x% isec3   ! (2)   ! bitmap
    y% rsec3 = x% rsec3   ! (2)
    !-----------------
    ! section 4 (data)
    !-----------------
    y% isec4 = x% isec4                ! binary data
    y% klenp = x% klenp                ! number of elements in array RSEC4
    call reallocate_data (y, y% klenp) ! decoded Real data
    if (associated (x% rsec4)) y% rsec4 (:y% klenp) = x% rsec4 (:y% klenp)
  END SUBROUTINE assign_grib1
!------------------------------------------------------------------------------
  subroutine bcast_grib1 (grib, source, comm)
  type (t_grib1) ,intent(inout)          :: grib   ! GRIB record to broadcast
  integer        ,intent(in)             :: source ! index of sending PE
  integer        ,intent(in)   ,optional :: comm   ! communicator

    integer :: c

    c = dace% comm; if (present(comm)) c = comm
    if (size_grib == 0) size_grib = size (transfer (grib,(/' '/)))

    call p_bcast_derivedtype (grib, size_grib, source, c)
    call p_bcast_ptr         (grib% kgrib,     source, c)
    call p_bcast_ptr         (grib% rsec4,     source, c)

  end subroutine bcast_grib1
!==============================================================================
  function local_dwd_grib1 (grib) result (local_dwd)
  type (t_grib1) ,intent(in) :: grib
  logical                    :: local_dwd

    local_dwd = (grib% isec1% center      == WMO0_DWD     .or.  &
                 grib% isec1% center      == WMO0_COSMO   .or.  &
                 grib% isec1% center      == WMO0_MSWISS  .or.  &
                 grib% isec1% center      == WMO0_COMET ) .and. &
                (grib% isec1% sub_center  == 0   .or.           &
                 grib% isec1% sub_center  == 255     )    .and. &
                 grib% isec1% local_flag  ==   1          .and. &
                (grib% isec1% local_ident == 152 .or.           &
                 grib% isec1% local_ident == 153 .or.           &
                 grib% isec1% local_ident >= 252 .and.          &
                 grib% isec1% local_ident <= 254)

  end function local_dwd_grib1
!------------------------------------------------------------------------------
  function local_ecmwf_grib1 (grib) result (local_ecmwf)
  type (t_grib1) ,intent(in) :: grib
  logical                    :: local_ecmwf

    local_ecmwf = grib% isec1% center      == WMO0_ECMWF .and. &
!                 grib% isec1% sub_center  == xxx .and. &
                  grib% isec1% local_flag  ==   1

  end function local_ecmwf_grib1
!------------------------------------------------------------------------------
  function bufr_size_enc (grib, verbose, nsec)
  !-------------------------------------------
  ! crosscheck the buffer size (array kgrib)
  !   from the encoded data (content of kgrib)
  !-------------------------------------------
  type (t_grib1) ,intent(in)            :: grib
  integer        ,intent(in)            :: verbose
  integer        ,intent(out) ,optional :: nsec (0:5)
  integer                               :: bufr_size_enc
    integer                :: i, n
    integer  , pointer     :: kg (:)
    integer  , pointer     :: pg (:)
#if !defined(__GFORTRAN__) && !defined(_CRAYFTN)
    character, allocatable :: cg (:)    ! Pointer is even slower on SX-6
#endif
    bufr_size_enc = 0
    if (verbose > 1) then
      print *
      print *,'bufr_size_enc:'
      print *
      print *,'blen        =', grib% blen
      print *,'size(kgrib) =', size(grib% kgrib)
      print *,'kleng       =', grib% kleng
      print *,'kword       =', grib% kword
      print *
    endif
    if (associated (grib% kgrib)) then
      allocate (kg (4*size(grib% kgrib)))
#if defined (__SX__) ||  defined(__ibm__)
      ! This FTRACE_REGION confuses SX's ftrace with MPI!
      FTRACE_BEGIN("bufr_size_enc:alt_code")
      !--------------------------------
      ! Workaround for slow transfer ()
      !--------------------------------
      if (little()) then
         do i = 1, size(grib% kgrib)
            n = 4*(i-1)
            kg(n+1) = ibits (grib% kgrib(i), 0, 8)
            kg(n+2) = ibits (grib% kgrib(i), 8, 8)
            kg(n+3) = ibits (grib% kgrib(i),16, 8)
            kg(n+4) = ibits (grib% kgrib(i),24, 8)
         end do
      else
         do i = 1, size(grib% kgrib)
            n = 4*(i-1)
            kg(n+1) = ibits (grib% kgrib(i),24, 8)
            kg(n+2) = ibits (grib% kgrib(i),16, 8)
            kg(n+3) = ibits (grib% kgrib(i), 8, 8)
            kg(n+4) = ibits (grib% kgrib(i), 0, 8)
         end do
      end if
      FTRACE_END  ("bufr_size_enc:alt_code")
#elif defined(__GFORTRAN__) || defined(_CRAYFTN)
      ! Avoid creation of explicit temporaries:
      kg = ICHAR (transfer (grib% kgrib, (/"a"/)))
#else
      allocate (cg (4*size(grib% kgrib)))
      FTRACE_BEGIN("bufr_size_enc:transfer")
      cg = transfer (grib% kgrib, cg)
      FTRACE_END  ("bufr_size_enc:transfer")

      FTRACE_BEGIN("bufr_size_enc:ichar")
      kg = ICHAR (cg)
      FTRACE_END  ("bufr_size_enc:ichar")
      deallocate (cg)
#endif
      !
      ! section 0
      !
      i = 0
      pg => kg
      if (size(pg) < 11) return
      if (verbose > 1) then
        write(6,'(a,i12,1x,a)')'s0,o5-7 :',256**2*pg(5)+256*pg(6)+pg(7),&
                                                              ' total length'
        write(6,'(a,i12,1x,a)')'s0,o8   :',pg(8),              ' edition'
        write(6,'(a,3i4,1x,a)')'s1,o1-3 :',pg(9),pg(10),pg(11),' PDS length'
      endif
      i = i + 8
      if (verbose > 0) then
        print *,'kgrib: section 0 =', 8,', sum =', i
      endif
      if (present(nsec)) nsec(0) = 8
      !
      ! section 1
      !
      pg => kg (i+1:)
      if (size(pg) < 41) return
      n = 256**2*pg(1)+256*pg(2)+pg(3)
      if (verbose > 1) then
        write(6,'(a,i12,1x,a)')'s1,o1-3  :', n       ,' PDS length'
        write(6,'(a,i12,1x,a)')'s1,o4    :',pg(4)    ,' table'
        write(6,'(a,i12,1x,a)')'s1,o5    :',pg(5)    ,' center'
        write(6,'(a,i12,1x,a)')'s1,o6    :',pg(6)    ,' process'
        write(6,'(a,i12,1x,a)')'s1,o7    :',pg(7)    ,' grid'
        write(6,'(a,i12,1x,a)')'s1,o8    :',pg(8)    ,' GDS/BDS present'
        write(6,'(a,i12,1x,a)')'s1,o9    :',pg(9)    ,' code'
        write(6,'(a,i12,1x,a)')'s1,o10   :',pg(10)   ,' leveltype'
        write(6,'(a,2i6,1x,a)')'s1,o11-12:',pg(11:12),' levels'
        write(6,'(a,5i2,1x,a)')'s1,o13-17:',pg(13:17),'   ccyymmddhh'
        write(6,'(a,i12,1x,a)')'s1,o18   :',pg(18)   ,' fc time unit'
        write(6,'(a,2i6,1x,a)')'s1,o19-20:',pg(19:20),' P1, P2'
        write(6,'(a,i12,1x,a)')'s1,o21   :',pg(21)   ,' time range'
        write(6,'(a,2i6,1x,a)')'s1,o22-23:',pg(22:23),' no. in average'
        write(6,'(a,i12,1x,a)')'s1,o24   :',pg(24)   ,' no. missing'
        write(6,'(a,i12,1x,a)')'s1,o25   :',pg(25)   ,' cc ref.time'
        write(6,'(a,i12,1x,a)')'s1,o26   :',pg(26)   ,' subcenter'
        write(6,'(a,2i6,1x,a)')'s1,o27-28:',pg(27:28),' scale factor'
        write(6,'(a,i12,1x,a)')'s1,o41   :',pg(41),' local extension identifier'
      endif
      if (iand(pg(8),64) /= 0)  &
       call finish('bufr_size_enc','section 3 present')
      if (iand(pg(8),128) == 0) &
       call finish('bufr_size_enc','section 2 not present')
      if (pg(8) /= 128)         &
       call finish('bufr_size_enc','section 2 not present or 3 present')
      i = i + n
      if (verbose > 0) then
        print *,'kgrib: section 1 =', n,', sum =', i
      endif
      if (present(nsec)) nsec(1) = n
      !
      ! section 2
      !
      pg => kg (i+1:)
      if (size(pg) < 6) return
      n = 256**2*pg(1)+256*pg(2)+pg(3)
      if (verbose > 1) then
        write(6,'(a,i12,1x,a)')'s2,o01-03:', n,' GDS length'
        write(6,'(a,i12,1x,a)')'s2,o04   :',pg(4),'NV'
        write(6,'(a,i12,1x,a)')'s2,o05   :',pg(5),'PV'
        write(6,'(a,i12,1x,a)')'s2,o06   :',pg(6),'data representation'
        !------------
        ! icosahedral
        !------------
        if (pg(6) == 192) then
        write(6,'(a,i12,1x,a)')'s2,o07-08:',pg(07)*256+pg(08),'fac2'
        write(6,'(a,i12,1x,a)')'s2,o09-10:',pg(09)*256+pg(10),'fac3'
        write(6,'(a,i12,1x,a)')'s2,o11-13:',pg(11)*256*256+&
                                            pg(12)*256+pg(13),'nd'
        write(6,'(a,i12,1x,a)')'s2,o14-16:',pg(14)*256*256+&
                                            pg(15)*256+pg(16),'ni'
        write(6,'(a,i12,1x,a)')'s2,o17   :',pg(17),'resolution,component flags'
        write(6,'(a,i12,1x,a)')'s2,o18-20:',pg(18)*256*256+&
                                            pg(19)*256+pg(20),'LatPol'
        write(6,'(a,i12,1x,a)')'s2,o21-23:',pg(21)*256*256+&
                                            pg(22)*256+pg(23),'LonPol'
        write(6,'(a,i12,1x,a)')'s2,o24-25:',pg(24)*256+pg(25),'LonDia'
        write(6,'(a,i12,1x,a)')'s2,o26-27:',pg(26)*256+pg(27),'notused'
        write(6,'(a,i12,1x,a)')'s2,o28   :',pg(28),'scanning mode'
        write(6,'(a,i12,1x,a)')'s2,o29   :',pg(29),'reserved'
        write(6,'(a,i12,1x,a)')'s2,o30   :',pg(30),'reserved'
        write(6,'(a,i12,1x,a)')'s2,o31   :',pg(31),'reserved'
        write(6,'(a,i12,1x,a)')'s2,o32   :',pg(32),'reserved'
        else
        !------------------------
        ! lat-lon incl. Gaussian:
        !------------------------
        write(6,'(a,i12,1x,a)')'s2,o07-08:',pg( 7)*256+pg( 8),'ni'
        write(6,'(a,i12,1x,a)')'s2,o09-10:',pg( 9)*256+pg(10),'nj'
        write(6,'(a,i12,1x,a)')'s2,o11-13:',pg(11)*256*256+&
                                            pg(12)*256+pg(13),'La1'
        write(6,'(a,i12,1x,a)')'s2,o14-16:',pg(14)*256*256+&
                                            pg(15)*256+pg(16),'Lo1'
        write(6,'(a,i12,1x,a)')'s2,o17   :',pg(17),'resolution,component flags'
        endif
      endif
      i = i + n
      if (verbose > 0) then
        print *,'kgrib: section 2 =', n,', sum =', i
      endif
      if (present(nsec)) nsec(2) = n
      !
      ! section 3
      !
      if (present(nsec)) nsec(3) = 0
      !
      ! section 4
      !
      pg => kg (i+1:)
      if (size(pg) < 14) return
      n = 256**2*pg(1)+256*pg(2)+pg(3)
      if (verbose > 1) then
        write(6,'(a,i12,1x,a)')'s4,o1-3  :', n,' BDS length'
        write(6,'(a,i12,1x,a)')'s4,o4    :',pg(4)/16,'flags'
        write(6,'(a,i12,1x,a)')'s4,o4    :',mod(pg(4),16),'unused bits'
        write(6,'(a,i12,1x,a)')'s4,o11   :',pg(11),'bits'
        write(6,'(a,i12,1x,a)')'s4,o14   :',pg(14)
      endif
      i = i + n
      if (verbose > 0) then
        print *,'kgrib: section 4 =', n,', sum =', i
      endif
      if (present(nsec)) nsec(4) = n
      !
      ! section 5
      !
      i = i + 4
      if (verbose > 0) then
        print *,'kgrib: section 5 =', 4,', sum =', i
      endif
      if (present(nsec)) nsec(5) = 4
      deallocate (kg)
      bufr_size_enc = i
    endif
  end function bufr_size_enc
!------------------------------------------------------------------------------
  function bufr_size_dec (grib, verbose, nsec)
  !-----------------------------------------------
  ! crosscheck the buffer size (array kgrib)
  !   from the decoded data (content of isec1,2,4)
  !-----------------------------------------------
  type (t_grib1) ,intent(in)            :: grib
  integer        ,intent(in)            :: verbose
  integer        ,intent(out) ,optional :: nsec(0:5)
  integer                               :: bufr_size_dec
    integer :: i, n
    logical :: reduced
    bufr_size_dec = -1
    !
    ! section 0
    !
    i = 8
    if (verbose > 1) then
      print *
      print *,'bufr_size_dec:'
      print *
    endif
    if (verbose > 0) then
      print *,'isec0: section 0 =', 8,', sum =', i
    endif
    if (present(nsec)) nsec(0) = 8
    !
    ! section 1
    !
    if (verbose > 1) then
      print *,'section 1'
      print *,'present_2_3    :', grib% isec1% present_2_3
      print *,'center         :', grib% isec1% center
      print *,'sub_center     :', grib% isec1% sub_center
      print *,'local_flag     :', grib% isec1% local_flag
      print *,'local_ident    :', grib% isec1% local_ident
      print *,'DWD extension  :', local_dwd  (grib)
      print *,'ECMWF extension:', local_ecmwf(grib)
    endif
    n = 28
    if (grib% isec1% local_flag /= 0) then
      select case (grib% isec1% local_ident)
      case (153) ! DWD local extension (Big-Ensemble forecasts)
        n = 70   ! 28 + 37 + 4 + 1 blank (padding) octet for even buffer size
      case (253) ! DWD local extension (Ensemble forecasts)
       select case (grib% isec1% center)
!      case (WMO0_MSWISS)
!       n = 64   ! 28 + 36     MeteoSwiss
       case default
        n = 66   ! 28 + 37 + 1 blank (padding) octet for even buffer size
       end select
      case (254) ! DWD local extension
        n = 54   ! 28 + 26
      case (1)   ! ECMWF local extension  1 - MARS labeling or ensemble forecst
        n = 52   !
        if (verbose > 1) then
#if defined (GRIBLIB_VERSION) && GRIBLIB_VERSION < 105
               print *, "Ensemble fcst. :", grib% isec1% local_use((51-42)+38)
#else
               print *, "Ensemble fcst. :", grib% isec1% local_use(43)
#endif
        end if
      case (14)  ! ECMWF local extension 14 - brightness temperatures
        n = 1080
      end select
    endif
    i = i + n
    if (verbose > 0) then
      print *,'isec1: section 1 =', n,', sum =', i
    endif
    if (present(nsec)) nsec(1) = n
    !
    ! section 2
    !
    select case (grib% isec2 (1))               ! Grid representation
    case default
       reduced = .false.
    case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
       reduced = (grib% isec4% n_data /= grib% isec4% non_miss)
       if (reduced .and. grib% isec4% non_miss == 0) then
          call message('bufr_size_dec','reduced grid, but isec4% non_miss==0?')
          call message('bufr_size_dec','treating as unreduced grid.')
          reduced = .false.
       end if
    end select

    if (verbose > 1) then
      print *,'section 2'
      print *,'repres    :',grib% isec2 (1)
      print *,'ni    ( 2):',grib% isec2 (2)     ! Number of longitudes
      print *,'nj    ( 3):',grib% isec2 (3)     ! Number of latitudes
      print *,'resol ( 6):',grib% isec2 (6)     ! Resolution & component flags
      select case (grib% isec2 (1))
      case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
         print *,'scan  (11):',grib% isec2 (11) ! Scanmode
      end select
      print *,'nv    (12):',grib% isec2 (12)
      print *,'reduc (17):',grib% isec2 (17)    ! =1: reduced grid
      if (reduced) then
         print *,'pl(1) (23):',grib% isec2 (23)     ! Longitudes near 1st pole
         print *,'pl_max    :',maxval (grib% isec2 (23:22+grib% isec2 (3)))
         print *, "(s4: non_miss /= n_data: assuming reduced grid)"
      end if
    endif
    n = grib% isec2 (12) * 4            ! List of vertical coordinate params.
    select case (grib% isec2 (1))
    case (WMO6_ROTLL)
       n = n + 42
    case default
       n = n + 32
    end select
    if (grib% isec2(17) == 1 .or. reduced) then
       n = n + 2 * grib% isec2(3)       ! Numbers of points/row of reduced grid
    end if
    i = i + n
    if (verbose > 0) then
      print *,'isec2: section 2 =', n,', sum =', i
    endif
    if (present(nsec)) nsec(2) = n
    !
    ! section 3
    !
    if (present(nsec)) nsec(3) = 0
    !
    ! section 4
    !
    if (verbose > 1) then
      print *,'section 4'
      print *,'n_data    :',grib% isec4% n_data
      print *,'bits      :',grib% isec4% bits
      print *,'grid_type :',grib% isec4% grid_type
      print *,'packing   :',grib% isec4% packing
      print *,'data_repr :',grib% isec4% data_repr
      print *,'flags     :',grib% isec4% flags
      print *,'reserved  :',grib% isec4% reserved
      print *,'matrix    :',grib% isec4% matrix
      print *,'bitmap2   :',grib% isec4% bitmap2
      print *,'width     :',grib% isec4% width
      print *,'bits2     :',grib% isec4% bits2
      print *,'wmo       :',grib% isec4% wmo
      print *,'start_cplx:',grib% isec4% start_cplx
      print *,'scale_cplx:',grib% isec4% scale_cplx
      print *,'j_cplx    :',grib% isec4% j_cplx
      print *,'k_cplx    :',grib% isec4% k_cplx
      print *,'m_cplx    :',grib% isec4% m_cplx
      print *,'non_miss  :',grib% isec4% non_miss
      print *,'reserved2 :',grib% isec4% reserved2
      print *,'o_coded   :',grib% isec4% o_coded
    endif
    if (reduced) then
       n = grib% isec4% non_miss * grib% isec4% bits / 8 + 11
    else
       n = grib% isec4% n_data   * grib% isec4% bits / 8 + 11
    end if
    if (mod (n, 2) == 1) n = n + 1          ! Pad to even size
    i = i + n
    if (verbose > 0) then
      print *,'isec4: section 4 =', n,', sum =', i
    endif
    if (present(nsec)) nsec(4) = n
    !
    ! section 5
    !
    if (verbose > 1) then
      print *,'section 4'
    endif
    n = 4
    i = i + n
    if (present(nsec)) nsec(5) = n
    if (verbose > 0) then
      print *,'isec : section 5 =', n,', sum =', i
    endif
    bufr_size_dec = i
  end function bufr_size_dec
!------------------------------------------------------------------------------
  subroutine get_octets (grib)
  type (t_grib1) ,intent(inout) :: grib

    integer :: i, j, n, l(0:5), o(grib% kword*4)

    j = 0
    if (little()) then
      do i = 1, grib% kword
        j = (i-1)*4
        o (j+1) = iand (      grib% kgrib(i)     , 255)
        o (j+2) = iand (ishft(grib% kgrib(i), -8), 255)
        o (j+3) = iand (ishft(grib% kgrib(i),-16), 255)
        o (j+4) = iand (ishft(grib% kgrib(i),-24), 255)
      end do
    else
      do i = 1, grib% kword
        j = (i-1)*4
        o (j+4) = iand (      grib% kgrib(i)     , 255)
        o (j+3) = iand (ishft(grib% kgrib(i), -8), 255)
        o (j+2) = iand (ishft(grib% kgrib(i),-16), 255)
        o (j+1) = iand (ishft(grib% kgrib(i),-24), 255)
      end do
    endif

    if (associated(grib% octets)) deallocate (grib% octets)
    allocate (grib% octets (grib% blen))
    grib% octets = o(1:grib% blen)

    n = bufr_size_enc (grib, 0, l)
    i = 0
    grib% os0 => grib% octets (i+1:i+l(0))
    i = i + l(0)
    grib% os1 => grib% octets (i+1:i+l(1))
    i = i + l(1)
    grib% os2 => grib% octets (i+1:i+l(2))
    i = i + l(2)
    grib% os3 => grib% octets (i+1:i+l(3))
    i = i + l(3)
    grib% os4 => grib% octets (i+1:i+l(4))
    i = i + l(4)
    grib% os5 => grib% octets (i+1:i+l(5))

  end subroutine get_octets
!------------------------------------------------------------------------------
  subroutine put_octets (grib)
  type (t_grib1) ,intent(inout) :: grib

    integer :: i, j, o (grib% kword*4)

    o(1:grib% blen) = grib% octets
    o ( grib% blen+1:) = 0

    j = 0
    if (little()) then
      do i = 1, grib% kword
        j = (i-1)*4
        grib% kgrib(i) =       o (j+1)     + &
                         ishft(o (j+2), 8) + &
                         ishft(o (j+3),16) + &
                         ishft(o (j+4),24)
      end do
    else
      do i = 1, grib% kword
        j = (i-1)*4
        grib% kgrib(i) =       o (j+4)     + &
                         ishft(o (j+3), 8) + &
                         ishft(o (j+2),16) + &
                         ishft(o (j+1),24)
      end do
    endif
    grib% kgrib( grib% kword+1:) = 0

  end subroutine put_octets
!------------------------------------------------------------------------------
  function diff_grib (grib1, grib2, verbose) result (diff)
  type (t_grib1) ,intent(inout) :: grib1
  type (t_grib1) ,intent(inout) :: grib2
  integer        ,intent(in)    :: verbose
  logical                       :: diff
  !======================================
  ! compares 2 GRIB records
  !
  ! verbose = 0: no output
  !           1: no decoding
  !           2: decode record
  !=======================================
    target  :: grib1,      grib2
    integer :: blen1,      blen2
    integer :: slen1(0:5), slen2(0:5)
    integer :: kleng
    integer :: blen
    integer :: n
    integer :: i, i1, i2, is
    type (t_sec0)      ,pointer :: s0_1, s0_2
    type (t_sec1)      ,pointer :: s1_1, s1_2
    type (t_s1_dwd)    ,pointer :: s1d1, s1d2
    type (t_s2_latlon) ,pointer :: s2l1, s2l2
    type (t_rsec2)     ,pointer :: s2r1, s2r2
    type (t_sec4)      ,pointer :: s4_1, s4_2
    character      ,allocatable :: c1(:), c2(:)

    diff = .false.

    !-----------------
    ! length of record
    !-----------------
    if (grib1% blen /= grib2% blen) then
      call header
      write (6,*) 'blen               :',grib1% blen, grib2% blen
    endif
    if (grib1% kleng /= grib2% kleng) then
      call header
      write (6,*) 'kleng              :',grib1% kleng, grib2% kleng
    endif
    if (grib1% kword /= grib2% kword) then
      call header
      write (6,*) 'kword              :',grib1% kword, grib2% kword
    endif
    if (diff .and. verbose == 0) return
    !--------------------------
    ! content of encoded record
    !--------------------------
    kleng = min (grib1% kleng, grib2% kleng)
    n = count (grib1% kgrib(:kleng) /= grib2% kgrib(:kleng))
    if (n > 0) then
      call header
      write (6,*) 'different words    :',n
      do i=1,kleng
        if (grib1% kgrib(i)/=grib2% kgrib(i)) then
          write (6,*) '1st diff in word   :',i,' of ',kleng,&
                      ', values=',grib1% kgrib(i),grib2% kgrib(i)
          exit
        end if
      end do
    endif
    if (.not. diff) return
    !------------------
    ! sizes of sections
    !------------------
    blen1 = bufr_size_enc (grib1, 0, slen1)
    blen2 = bufr_size_enc (grib2, 0, slen2)
    if (blen1/=blen2 .or. any (slen1/=slen2)) then
      call header
      write (6,*) 'sum(len(secs))     :',blen1,    blen2
      write (6,*) 'len(sec 0)         :',slen1(0), slen2(0)
      write (6,*) 'len(sec 1)         :',slen1(1), slen2(1)
      write (6,*) 'len(sec 2)         :',slen1(2), slen2(2)
      write (6,*) 'len(sec 3)         :',slen1(3), slen2(3)
      write (6,*) 'len(sec 4)         :',slen1(4), slen2(4)
      write (6,*) 'len(sec 5)         :',slen1(5), slen2(5)
    endif
    !-------------------
    ! octets in sections
    !-------------------
    if (diff) then
      allocate (c1(max(grib1% kleng, grib2% kleng)*4))
      allocate (c2(max(grib1% kleng, grib2% kleng)*4))
      c1 = transfer (grib1% kgrib, c1)
      c2 = transfer (grib2% kgrib, c2)
      i1 = 0; i2 = 0
      do i=0,5
        blen = min (slen1(i),slen2(i))
        n = count (c1(i1+1:i1+blen) /= c2(i2+1:i2+blen))
        if (n>0) then
          write (6,'(a,i6,a,i2)') 'different octets   :',n,' in section ',i
          do is=1,blen
            if (c1(i1+is) /= c2(i2+is)) then
              write (6,'(a,i6,a,i6,a,2i4)')          &
               '1st diff in oct.   :',is,' of ',blen,&
               ', values=',ICHAR(c1(i1+is)),         &
                           ICHAR(c2(i2+is))
!             exit
            end if
          end do
        endif

        n = slen1(i) - slen2(i)
        if (n>0) then
          write (6,'(a,i6,a,i2)') 'additional octets  :',n,' in section ',i
          do is=blen+1, blen+n
            write (6,'(a,i6,a,i6,a,i4,a)')            &
              'additional  oct.   :',is,' of ',blen+n,&
              ', values=',ICHAR(c1(i1+is)),' ---'
          end do
        end if

        n = slen2(i) - slen1(i)
        if (n>0) then
          write (6,'(a,i6,a,i2)') 'additional octets  :',n,' in section ',i
          do is=blen+1, blen+n
            write (6,'(a,i6,a,i6,a,a,i4)')            &
              'additional  oct.   :',is,' of ',blen+n,&
              ', values=',' ---',ICHAR(c2(i2+is))
          end do
        end if

        i1 = i1 + slen1(i)
        i2 = i2 + slen2(i)
      end do
      deallocate (c1, c2)
    endif
    !------------
    ! decode data
    !------------
    if (verbose <= 1) return
    call gribex_90 (grib1,'D')
    call gribex_90 (grib2,'D')
    !----------
    ! section 0
    !----------
    write(6,*)
    write(6,*) 'Diffs in Section 0 (Indicator Section):'
    write(6,*)
    s0_1 => grib1% isec0
    s0_2 => grib2% isec0
    call diffi (s0_1% n_octets  , s0_2% n_octets  , 'n_octets           :')
    call diffi (s0_1% edition   , s0_2% edition   , 'edition            :')
    !----------
    ! section 1
    !----------
    write(6,*)
    write(6,*) 'Diffs in Section 1 (Product Definition Section):'
    write(6,*)
    s1_1 => grib1% isec1
    s1_2 => grib2% isec1
    call diffi (s1_1% table      , s1_2% table      , 'table              :')
    call diffi (s1_1% center     , s1_2% center     , 'center             :')
    call diffi (s1_1% process    , s1_2% process    , 'process            :')
    call diffi (s1_1% grid       , s1_2% grid       , 'grid               :')
    call diffi (s1_1% present_2_3, s1_2% present_2_3, 'present_2_3        :')
    call diffi (s1_1% code       , s1_2% code       , 'code               :')
    call diffi (s1_1% level_type , s1_2% level_type , 'level_type         :')
    call diffi (s1_1% level_st   , s1_2% level_st   , 'level_st           :')
    call diffi (s1_1% level_b    , s1_2% level_b    , 'level_b            :')
    call diffi (s1_1% year       , s1_2% year       , 'year               :')
    call diffi (s1_1% month      , s1_2% month      , 'month              :')
    call diffi (s1_1% day        , s1_2% day        , 'day                :')
    call diffi (s1_1% hour       , s1_2% hour       , 'hour               :')
    call diffi (s1_1% minute     , s1_2% minute     , 'minute             :')
    call diffi (s1_1% time_unit  , s1_2% time_unit  , 'time_unit          :')
    call diffi (s1_1% p1         , s1_2% p1         , 'p1                 :')
    call diffi (s1_1% p2         , s1_2% p2         , 'p2                 :')
    call diffi (s1_1% time_range , s1_2% time_range , 'time_range         :')
    call diffi (s1_1% n_average  , s1_2% n_average  , 'n_average          :')
    call diffi (s1_1% n_missing  , s1_2% n_missing  , 'n_missing          :')
    call diffi (s1_1% century    , s1_2% century    , 'century            :')
    call diffi (s1_1% sub_center , s1_2% sub_center , 'sub_center         :')
    call diffi (s1_1% factor     , s1_2% factor     , 'factor             :')
    call diffi (s1_1% local_flag , s1_2% local_flag , 'local_flag         :')
    call diffii(s1_1% reserved   , s1_2% reserved   , 'reserved           :' ,24)
    call diffi (s1_1% local_ident, s1_2% local_ident, 'local_ident        :')
    call diffii(s1_1% local_use  , s1_2% local_use  , 'local_use          :' ,37)
    !-------------------------
    ! section 1, DWD extension
    !-------------------------
    if (local_dwd (grib1) .and. local_dwd (grib2)) then
    write(6,*)
    write(6,*) 'Diffs in Section 1 (DWD Extension):'
    write(6,*)
    s1d1 => grib1% s1_dwd
    s1d2 => grib2% s1_dwd
    call diffi (s1d1% local_ident  , s1d2% local_ident  , 'local_ident        :')
    call diffi (s1d1% day_number   , s1d2% day_number   , 'day_number         :')
    call diffi (s1d1% record_number, s1d2% record_number, 'record_number      :')
    call diffi (s1d1% decoding     , s1d2% decoding     , 'decoding           :')
    call diffi (s1d1% element_no   , s1d2% element_no   , 'element_no         :')
    call diffi (s1d1% year         , s1d2% year         , 'year               :')
    call diffi (s1d1% month        , s1d2% month        , 'month              :')
    call diffi (s1d1% day          , s1d2% day          , 'day                :')
    call diffi (s1d1% hour         , s1d2% hour         , 'hour               :')
    call diffi (s1d1% minute       , s1d2% minute       , 'minute             :')
    call diffi (s1d1% exp          , s1d2% exp          , 'exp                :')
    call diffi (s1d1% run_type     , s1d2% run_type     , 'run_type           :')
    call diffi (s1d1% user_id      , s1d2% user_id      , 'user_id            :')
    call diffi (s1d1% experiment_id, s1d2% experiment_id, 'experiment_id      :')
    call diffi (s1d1% ensemble_id  , s1d2% ensemble_id  , 'ensemble_id        :')
    call diffi (s1d1% ensemble_size, s1d2% ensemble_size, 'ensemble_size      :')
    call diffi (s1d1% ensemble_no  , s1d2% ensemble_no  , 'ensemble_no        :')
    call diffi (s1d1% major_version, s1d2% major_version, 'major_version      :')
    call diffi (s1d1% minor_version, s1d2% minor_version, 'minor_version      :')
    endif
    !----------
    ! section 2
    !----------
    write(6,*)
    write(6,*) 'Diffs in Section 2 (Grid Definition Section):'
    write(6,*)
    call diffi (grib1% isec2(1)    , grib2% isec2(1)    , 'representation     :')
    call diffii(grib1% isec2(2:)   , grib2% isec2(2:)   , 'remaining section 2:' ,1)
    if (grib1% isec2(1) == grib2% isec2(1)) then
    select case (grib1% isec2(1))
    !-------------------
    ! section 2 (latlon)
    !-------------------
    case (WMO6_LATLON, WMO6_ROTLL)
    write(6,*)
    write(6,*) '  (rotated) lat-lon grid :'
    write(6,*)
    s2l1 => grib1% latlon
    s2l2 => grib2% latlon
    call diffi (s2l1% repr         , s2l2% repr         , 'repr               :')
    call diffi (s2l1% ni           , s2l2% ni           , 'ni                 :')
    call diffi (s2l1% nj           , s2l2% nj           , 'nj                 :')
    call diffi (s2l1% lat_first    , s2l2% lat_first    , 'lat_first          :')
    call diffi (s2l1% lon_first    , s2l2% lon_first    , 'lon_first          :')
    call diffi (s2l1% increments   , s2l2% increments   , 'increments         :')
    call diffi (s2l1% lat_last     , s2l2% lat_last     , 'lat_last           :')
    call diffi (s2l1% lon_last     , s2l2% lon_last     , 'lon_last           :')
    call diffi (s2l1% di           , s2l2% di           , 'di                 :')
    call diffi (s2l1% dj           , s2l2% dj           , 'dj                 :')
    call diffi (s2l1% scan_mode    , s2l2% scan_mode    , 'scan_mode          :')
    call diffi (s2l1% nvcp         , s2l2% nvcp         , 'nvcp               :')
    call diffi (s2l1% lat_rot      , s2l2% lat_rot      , 'lat_rot            :')
    call diffi (s2l1% lon_rot      , s2l2% lon_rot      , 'lon_rot            :')
    call diffi (s2l1% lat_strech   , s2l2% lat_strech   , 'lat_strech         :')
    call diffi (s2l1% lon_strech   , s2l2% lon_strech   , 'lon_strech         :')
    call diffi (s2l1% reduced      , s2l2% reduced      , 'reduced            :')
    call diffi (s2l1% earth        , s2l2% earth        , 'earth              :')
    call diffi (s2l1% components   , s2l2% components   , 'components         :')
    call diffii(s2l1% reserved     , s2l2% reserved     , 'reserved           :' ,19)
!   INTEGER            :: isec2  (1024) ! Integer data
!   TYPE (t_rsec2)     :: rsec2 ! (512) ! Real data
!   TYPE (t_s2_latlon) :: latlon        ! Integer data, lat/lon    grid
!   TYPE (t_s2_gauss)  :: gauss         ! Integer data, Gauss      grid
!   TYPE (t_s2_tri)    :: tri           ! Integer data, triangular grid
!   TYPE (t_s2_sph)    :: sph           ! Integer data, sph. harm. grid
    case default
    write(6,*)
    write(6,*) '  unkown representation:', grib1% isec2(1)
    write(6,*)
    end select
    endif

    !---------------------
    ! section 2, real part
    !---------------------
    write(6,*)
    write(6,*) '  section 2, real part :'
    write(6,*)
    s2r1 => grib1% rsec2
    s2r2 => grib2% rsec2
    call diffr (s2r1% rot_angle , s2r2% rot_angle , 'rot_angle          :')
    call diffr (s2r1% str_factor, s2r2% str_factor, 'str_factor         :')
    call diffrr(s2r1% reserved  , s2r2% reserved  , 'reserved           :',  0)
    call diffrr(s2r1% vcp       , s2r2% vcp       , 'vcp                :',  0)

    !------------------------
    ! section 4, integer part
    !------------------------
    write(6,*)
    write(6,*) 'Diffs in Section 4 (Binary Data Section):'
    write(6,*)
    s4_1 => grib1% isec4
    s4_2 => grib2% isec4
    call diffi (s4_1% n_data    , s4_2% n_data    , 'n_data             :')
    call diffi (s4_1% bits      , s4_2% bits      , 'bits               :')
    call diffi (s4_1% grid_type , s4_2% grid_type , 'grid_type          :')
    call diffi (s4_1% packing   , s4_2% packing   , 'packing            :')
    call diffi (s4_1% data_repr , s4_2% data_repr , 'data_repr          :')
    call diffi (s4_1% flags     , s4_2% flags     , 'flags              :')
    call diffi (s4_1% reserved  , s4_2% reserved  , 'reserved           :')
    call diffi (s4_1% matrix    , s4_2% matrix    , 'matrix             :')
    call diffi (s4_1% bitmap2   , s4_2% bitmap2   , 'bitmap2            :')
    call diffi (s4_1% width     , s4_2% width     , 'width              :')
    call diffi (s4_1% bits2     , s4_2% bits2     , 'bits2              :')
    call diffii(s4_1% wmo       , s4_2% wmo       , 'wmo                :', 11) !12:15
    call diffi (s4_1% start_cplx, s4_2% start_cplx, 'start_cplx         :')
    call diffi (s4_1% scale_cplx, s4_2% scale_cplx, 'scale_cplx         :')
    call diffi (s4_1% j_cplx    , s4_2% j_cplx    , 'j_cplx             :')
    call diffi (s4_1% k_cplx    , s4_2% k_cplx    , 'k_cplx             :')
    call diffi (s4_1% m_cplx    , s4_2% m_cplx    , 'm_cplx             :')
    call diffi (s4_1% non_miss  , s4_2% non_miss  , 'non_miss           :')
    call diffii(s4_1% reserved2 , s4_2% reserved2 , 'reserved2          :', 21) !22:33
    call diffi (s4_1% o_coded   , s4_2% o_coded   , 'o_coded            :')
    call diffii(s4_1% remaining , s4_2% remaining , 'remaining          :', 34) !35:512
    !---------------------
    ! section 4, real part
    !---------------------
    call diffrr(grib1% rsec4    , grib2% rsec4    , 'rsec4              :',  0)

print *,'xxxxxxxxxxxxxxxxxxxxx'

  contains

    subroutine header
      if (.not. diff .and. verbose > 0) then
        write (6,*)
        write (6,*) 'GRIB buffers differ:'
        write (6,*) 'krec               :',grib1% krec, grib2% krec
        write (6,*) 'kpos               :',grib1% kpos, grib2% kpos
      endif
      diff = .true.
    end subroutine header

    subroutine diffi (i1, i2, name)
    integer          ,intent(in) :: i1, i2
    character(len=*) ,intent(in) :: name
      if (i1 /= i2) write (6,*) name, i1,i2
    end subroutine diffi

    subroutine diffii (i1, i2, name, io)
    integer          ,intent(in) :: i1(:), i2(:), io
    character(len=*) ,intent(in) :: name
      integer :: n, m, i
      m = min (size(i1), size(i2))
      n = count (i1(1:m)/=i2(:m))
      if (n > 0) then
!
! inlining for compiler error in icon code:
! ftn-2116 crayftn: INTERNAL
!   "/opt/cray/cce/8.4.1/cftn/x86-64/lib/optcg" was terminated due to receipt of signal 013:  Segmentation fault.
!
!       call header

        if (.not. diff .and. verbose > 0) then
          write (6,*)
          write (6,*) 'GRIB buffers differ:'
          write (6,*) 'krec               :',grib1% krec, grib2% krec
          write (6,*) 'kpos               :',grib1% kpos, grib2% kpos
        endif
        diff = .true.

        write (6,*) name, n,' different words of ',m
        do i=1,m
          if (i1(i)/=i2(i)) then
            write (6,*) '1st diff in word   :',i+io,&
                        ', values=',i1(i),i2(i)
            exit
          end if
        end do
      endif
    end subroutine diffii

    subroutine diffr (r1, r2, name)
    real(wp)         ,intent(in) :: r1, r2
    character(len=*) ,intent(in) :: name
      if (r1 /= r2) write (6,*) name, r1, r2
    end subroutine diffr

    subroutine diffrr (r1, r2, name, io)
    real(wp)         ,intent(in) :: r1(:), r2(:)
    integer          ,intent(in) :: io
    character(len=*) ,intent(in) :: name
      integer :: n, m, i
      m = min (size(r1), size(r2))
      n = count (r1(1:m)/=r2(:m))
      if (n > 0) then
        write (6,*) name, n,' different words of ',m
        do i=1,m
          if (r1(i)/=r2(i)) then
            write (6,*) '1st diff in word   :',i+io,&
                        ', values=',r1(i),r2(i)
            exit
          end if
        end do
      endif
    end subroutine diffrr

  end function diff_grib
!------------------------------------------------------------------------------
  subroutine print_isec4 (s4i)
  type (t_sec4) ,intent(in) :: s4i
  !--------------------------------
  ! print integer data of section 4
  !--------------------------------
    write (6,'(a)')
    write (6,'(a)') ' Section 1 - Data Section.'
    write (6,'(a)') ' -------------------------'
    write (6,'(a)')
    write (6,'(a, i5,a)') ' n_data    :', s4i% n_data    ,' Number of data values in array PSEC4.'
    write (6,'(a, i5,a)') ' bits      :', s4i% bits      ,' Number of bits used for each encoded value.'
    write (6,'(a, i5,a)') ' grid_type :', s4i% grid_type ,' 0: Grid point data, 128: Sph.harm.coefficients'
    write (6,'(a, i5,a)') ' packing   :', s4i% packing   ,' 0: Simple packing,   64: Complex packing.'
    write (6,'(a, i5,a)') ' data_repr :', s4i% data_repr ,' 0: Floating point,   32: Integer data.'
    write (6,'(a, i5,a)') ' flags     :', s4i% flags     ,' 0: no                16: Additional flags.'
    write (6,'(a, i5,a)') ' reserved  :', s4i% reserved  ,' Reserved. Set to 0.'
    write (6,'(a, i5,a)') ' matrix    :', s4i% matrix    ,' 0: single datum      64: Matrix at each grid point.'
    write (6,'(a, i5,a)') ' bitmap2   :', s4i% bitmap2   ,' 0: no                32: Secondary bitmaps present.'
    write (6,'(a, i5,a)') ' width     :', s4i% width     ,' 0: constant width    16: different width of 2nd order values.'
    write (6,'(a, i5,a)') ' bits2     :', s4i% bits2     ,' Number of bits for second order values.'
    write (6,'(a,4i5,a)') ' wmo       :', s4i% wmo
    write (6,'(a, i5,a)') ' start_cplx:', s4i% start_cplx,' For complex packing, start of packed data values.'
    write (6,'(a, i5,a)') ' scale_cplx:', s4i% scale_cplx,' For complex packing, the scaling factor P.'
    write (6,'(a, i5,a)') ' j_cplx    :', s4i% j_cplx    ,' For complex packing, pentagonal resolution parameter J.'
    write (6,'(a, i5,a)') ' k_cplx    :', s4i% k_cplx    ,' For complex packing, pentagonal resolution parameter K.'
    write (6,'(a, i5,a)') ' m_cplx    :', s4i% m_cplx    ,' For complex packing, pentagonal resolution parameter M.'
    write (6,'(a, i5,a)') ' non_miss  :', s4i% non_miss  ,' The number of non-missing values.'
    write (6,'(a,2i5,a)') ' reserved2 :', minval(s4i% reserved2) ,maxval(s4i% reserved2)
    write (6,'(a, i5,a)') ' o_coded   :', s4i% o_coded   ,' offset bit pointer to coded values in GRIB record.'
    write (6,'(a)')
  end subroutine print_isec4
!==============================================================================
  elemental function jchar (c)
  integer               :: jchar
  character ,intent(in) :: c
  !------------------------------------------------------------------------
  ! On IBM the ichar intrinsic returns negative numbers.
  ! (The treatment of character codes > 127 is not defined by the standard)
  ! Thus we provide this replacement.
  !------------------------------------------------------------------------
    jchar = ichar (c)
    if (jchar < 0) jchar = jchar + 256
  end function jchar
!==============================================================================
END MODULE mo_emos_grib1
