!
!+ Handle ECMWF local GRIB extensions (MPIfB GRIB library)
!
! $Id$
!
MODULE mo_ecmwf_grib
!
! Description:
! ECMWF local GRIB extensions
! used by module mo_emos_grib1.f90  (MPIfM EMOS compatible GRIB library)
!
! Up to now only: Definition 14 - Brightness temperature
!                 Definition 1  - Mars labelling or ensemble forecast data
!                                 (partially implemented)
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004-2007  original code
! Harald Anlauf   DWD  2007       adjustments for GRIB library 1.0.7
!------------------------------------------------------------------------------
  !=============
  ! Modules used
  !=============
  USE mo_kind,       ONLY: i4      ! kind parameters
  USE mo_exception,  ONLY: finish  ! abort on error condition
  USE mo_endian,     ONLY: little  ! returns .true. for little endian
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: t_s1_ecmwf     ! ECMWF local extension data type
  public :: decode_ecmwf   ! decode the ECMWF local extensions
  public :: encode_ecmwf   ! encode the ECMWF local extensions
  public :: grprs1ec       ! print ECMWF local use part of section 1 (PDS)

  !---------------------------------------------------------------------
  ! Data type definition for GRIB PDS (section 1) ECMWF local extensions
  !---------------------------------------------------------------------
  type t_s1_ecmwf                   ! octets|isec1|
    integer          :: local_ident ! 41    | 37  | local GRIB use identifier.
    integer          :: class       ! 42    | 38  |
    integer          :: type        ! 43    | 39  |
    integer          :: stream      ! 44-45 | 40  |
    character(len=4) :: expid       ! 46-49 | 41  | version number or exp. id.
    integer          :: i_en        ! 50    | 42  | ensemble forecast number
    integer          :: n_en        ! 51    | 43  | total # of fc in ensemble
    !------------------------------------------------
    ! Definition 14 specific: Brightness temperatures
    !------------------------------------------------
    integer          :: channel     ! 52    | 44  | channel number
    integer          :: scaling     ! 53-56 | 45  | frequencies scaling factor
    integer          :: nf          ! 57    | 46  | total number of frequencies
    !                               | 58-60 |     | spare
!   integer          :: freq(255)   !61-1080| 47- | list of scaled frequencies
  end type t_s1_ecmwf

contains
!==============================================================================
  subroutine decode_ecmwf (ecmwf, ident, lu)
  !----------------------------------
  ! decode the ECMWF local extensions
  !----------------------------------
  type (t_s1_ecmwf) ,intent(out) :: ecmwf    ! ECMWF local extension data type
  integer           ,intent(in)  :: ident    ! local use definition identifier
  integer           ,intent(in)  :: lu (42:) ! octets
    character(len=4)   :: cext
    integer(i4)        :: id
    ecmwf% local_ident = ident
    ecmwf% class       = lu (42)
    ecmwf% type        = lu (43)
#if defined (GRIBLIB_VERSION) && GRIBLIB_VERSION < 105
    ! GRIB library versions < 1.0.5:
    ecmwf% stream      = lu (44)*256 + lu (45)
    if (little()) then
      id = ishft(lu(49),24) + ishft(lu(48),16) + ishft(lu(47),8) + lu(46)
    else
      id = ishft(lu(46),24) + ishft(lu(47),16) + ishft(lu(48),8) + lu(49)
    endif
    ecmwf% expid = transfer (id, ecmwf% expid)
    ecmwf% i_en        = lu (50)
    ecmwf% n_en        = lu (51)
#else
    ! GRIB library versions >= 1.0.5 decode the ECMWF extension 1 differently
    ecmwf% stream      = lu (44)
    id                 = int (lu(45), kind=i4)
    ecmwf% expid       = transfer (id, ecmwf% expid)    ! GRIB 1.0.7!
    ecmwf% i_en        = lu (46)
    ecmwf% n_en        = lu (47)
#endif
    !------------------------
    ! Definition 14 specific:
    !------------------------
    ecmwf% channel     = 0
    ecmwf% scaling     = 0
    ecmwf% nf          = 0
    select case (ident)
    case ( 1)
      !----------------------------------------
      ! MARS labeling or ensemble forecast data
      !----------------------------------------
    case (14)
      !-----------------------
      ! Brightness temperature
      !-----------------------
      call finish ("decode_ecmwf", &
           "local extension 14 not adapted to GRIB library >= 1.0.5")
      ecmwf% channel   =   lu(52)
      ecmwf% scaling   = ((lu(53)*256+lu(54))*256+lu(55))*256+lu(56)
      ecmwf% nf        =   lu(57)
    case default
      write (cext,'(i4)') ident
      call finish ('decode_ecmwf','local extension not implemented:'//cext)
    end select
  end subroutine decode_ecmwf
!------------------------------------------------------------------------------
  subroutine encode_ecmwf (ecmwf, lu)
  !----------------------------------
  ! decode the ECMWF local extensions
  !----------------------------------
  type (t_s1_ecmwf) ,intent(in)   :: ecmwf    ! ECMWF local extension data type
  integer           ,intent(out)  :: lu (42:) ! octets
    integer                       :: ident    ! local use definition identifier
    integer                       :: scaling  ! temporary
    character(len=4)              :: cext
    integer(i4)                   :: id
    ident   =          ecmwf% local_ident
    lu (42) =          ecmwf% class
    lu (43) =          ecmwf% type
#if defined (GRIBLIB_VERSION) && GRIBLIB_VERSION < 105
    ! GRIB library versions < 1.0.5:
    lu (44) =          ecmwf% stream / 256
    lu (45) =      mod(ecmwf% stream , 256)
    id = transfer (ecmwf% expid, id)
    if (little()) then
      lu(49) = iand (ishft(id,-24),255)
      lu(48) = iand (ishft(id,-16),255)
      lu(47) = iand (ishft(id,- 8),255)
      lu(46) = iand (      id     ,255)
    else
      lu(46) = iand (ishft(id,-24),255)
      lu(47) = iand (ishft(id,-16),255)
      lu(48) = iand (ishft(id,- 8),255)
      lu(49) = iand (      id     ,255)
    endif
    lu (50) =          ecmwf% i_en
    lu (51) =          ecmwf% n_en
    lu (52) =                 0
#else
    ! GRIB library versions >= 1.0.5 encode the ECMWF extension 1 differently
    lu (44) =          ecmwf% stream
    id = transfer (ecmwf% expid, id)
    lu (45) =                 id                        ! GRIB 1.0.7!
    lu (46) =          ecmwf% i_en
    lu (47) =          ecmwf% n_en
    lu (48) =                 0
#endif
    !------------------------
    ! Definition 14 specific:
    !------------------------
    select case (ident)
    case ( 1)
      !----------------------------------------
      ! MARS labeling or ensemble forecast data
      !----------------------------------------
    case (14)
      !-----------------------
      ! Brightness temperature
      !-----------------------
      call finish ("encode_ecmwf", &
           "local extension 14 not adapted to GRIB library >= 1.0.5")
      scaling = ecmwf% scaling
      lu(57)  = ecmwf% nf
      lu(56)  = mod(scaling , 256); scaling = ecmwf% scaling / 256
      lu(55)  = mod(scaling , 256); scaling = ecmwf% scaling / 256
      lu(54)  = mod(scaling , 256); scaling = ecmwf% scaling / 256
      lu(53)  =     scaling
      lu(52)  = ecmwf% channel
    case default
      write (cext,'(i4)') ident
      call finish ('encode_ecmwf','local extension not implemented:'//cext)
    end select
  end subroutine encode_ecmwf
!------------------------------------------------------------------------------
  subroutine grprs1ec (ecmwf)
  TYPE (t_s1_ecmwf) ,INTENT(in) :: ecmwf
  !----------------------------------------------
  ! print ECMWF local use part of section 1 (PDS)
  !----------------------------------------------
    if (ecmwf% local_ident > 0) then
      write(6,"()")
      write(6,"(a)")    ' Section 1 - Product Definition Section, Local Use.'
      write(6,"(a)")    ' --------------------------------------------------'
      write(6,"()")
      write(6,"(a39,i8)")'ECMWF local use identifier    ', ecmwf% local_ident
      write(6,"(a39,i8)")'class                         ', ecmwf% class
      write(6,"(a39,i8)")'type                          ', ecmwf% type
      write(6,"(a39,i8)")'stream                        ', ecmwf% stream
      write(6,"(a39,a8)")'version no or exp.id.         ', ecmwf% expid
      write(6,"(a39,i8)")'ensemble forecast number      ', ecmwf% i_en
      write(6,"(a39,i8)")'total number of fc in ensemble', ecmwf% n_en
      select case (ecmwf% local_ident)
      case ( 1)
        write(6,"(a40)")   'MARS labeling or ensemble forecast data!'
      case (14)
        write(6,"(a39)")   'Brightness temperature :'
        write(6,"(a39,i8)")'channel number               ', ecmwf% channel
        write(6,"(a39,i8)")'scaling factor               ', ecmwf% scaling
        write(6,"(a39,i8)")'total number of frequencies  ', ecmwf% nf
      case default
        write(6,"(a39)")   'unknown ECMWF local extension!'
      end select
    endif
  end subroutine grprs1ec
!==============================================================================
end module mo_ecmwf_grib
