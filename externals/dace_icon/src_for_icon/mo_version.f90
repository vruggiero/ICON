!
!+ Provide 3D-Var version number and modification date
!
MODULE mo_version
!
! Description:
!   Provide 3D-Var version number and modification date
!   (internal, not DWD SCCS)
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
! V1_2         2008/12/04 Andreas Rhodin
!  bug fix for previous release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_40        2015-02-27 Harald Anlauf
!  Diagnose compiler vendor and version (H.Anlauf)
! V1_42        2015-06-08 Harald Anlauf
!  Disentangle modules info_3dvar and mo_version
! V1_48        2016-10-06 Harald Anlauf
!  compiler_info: fix for obsolete Intel compilers (v13)
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007-2008
!------------------------------------------------------------------------------
!
! Old outdated labeling:
!
! Revision 1.7  2008/09/25 13:21:49  m214030
! remove trailing blanks
!
! Revision 1.6  2008/09/12 06:59:01  m214030
! add file headers for DWD SCCS
!
! Revision 1.5  2007/06/01 18:15:32  m214030
! version for ECMWF / Vietnam
!
! Revision 1.4  2007/04/26 08:33:03  m214030
! wavelet representation of horizontal B matrix
!
! Revision 1.3  2007/02/08 18:00:06  m214030
! use ECMWF fields, GPSRO operator 1,2,3, AIREP moisture
!
! Revision 1.2  2006/11/21 15:19:30  m214030
! 1D-GPS RO operator implemented
!
! Revision 1.1  2006/11/07 08:43:29  m214030
! first revision
!
!------------------------------------------------------------
#if defined (__FLANG) || defined (__NEC__) || defined (__INTEL_COMPILER)
  use, intrinsic :: iso_fortran_env, only: compiler_version
#endif

!use info_3dvar, only: info_getvalue
implicit none

private
public :: var3d_version  ! return version number         'X.YY'
public :: var3d_date     !        version date           'YYYY-MM-DD'
public :: var3d_time     !        version time           'HH:MM:SS'
public :: major          !        major version number     X
public :: minor          !        minor version number    YY
public :: compiler_info  ! return compiler vendor and version information

#include "info_3dvar.incf"

contains

!------------------------------------------------------------------------------
  function var3d_version ()
  character(len=8) :: var3d_version
  !-----------------------------
  ! return version number 'X.YY'
  !-----------------------------
    character(len=8) :: tmp
    integer          :: i

!   tmp = info_getvalue ('t')
    tmp = INFO_Tag
    i   = index (tmp,'_')
    tmp(i:i) = '.'
    var3d_version = tmp(2:)

  end function var3d_version
!------------------------------------------------------------------------------
  function var3d_date ()
  character(len=10) :: var3d_date
  !---------------------------------
  ! return version date 'YYYY-MM-DD'
  !---------------------------------
    character(len=19) :: tmp            ! 'YYYY-MM-DD HH:MM:SS'

!   tmp        = info_getvalue ('i')
    tmp        = INFO_Date
    var3d_date = tmp (1:10)

  end function var3d_date
!------------------------------------------------------------------------------
  function var3d_time ()
  character(len=8) :: var3d_time
  !-------------------------------
  ! return version time 'HH:MM:SS'
  !-------------------------------
    character(len=19) :: tmp            ! 'YYYY-MM-DD HH:MM:SS'

!   tmp        = info_getvalue ('i')
    tmp        = INFO_Date
    var3d_time = tmp (12:19)

  end function var3d_time
!------------------------------------------------------------------------------
  function major()
  integer :: major
  !----------------------------
  ! return major version number
  !----------------------------
    character(len=8) :: tmp, c
    integer          :: i
    integer, save    :: major_ = -1

    if (major_ < 0) then
!      tmp = info_getvalue ('t')
       tmp = INFO_Tag
       i   = index (tmp,'_')
       c   = tmp   (2:i-1)
       read (c,*) major_
    end if
    major = major_

  end function major
!------------------------------------------------------------------------------
  function minor()
  integer :: minor
  !----------------------------
  ! return minor version number
  !----------------------------
    character(len=8) :: tmp, c
    integer          :: i
    integer, save    :: minor_ = -1

    if (minor_ < 0) then
!      tmp = info_getvalue('t')
       tmp = INFO_Tag
       i   = index(tmp,'_')
       c   = tmp(i+1:)
       read (c,*) minor_
    end if
    minor = minor_

  end function minor
!------------------------------------------------------------------------------
  subroutine compiler_info (vendor, version, major, minor, patchlevel)
    character(*),intent(out), optional :: vendor, version
    integer     ,intent(out), optional :: major, minor, patchlevel
    !------------------------------------------------
    ! return compiler vendor, version/release string,
    ! major version, minor and patchlevel.
    !------------------------------------------------
    character(len=64)  :: vendor_
    character(len=140) :: version_
    integer            :: major_, minor_, patchlevel_, ios
#if defined (NAGFOR)
    integer            :: i
#endif
#if defined (__SUNPRO_F95) || defined (_CRAYFTN) || defined (__FLANG) || \
    defined (__NEC__) || defined (__INTEL_COMPILER)
    integer            :: i, j
    character(len=64)  :: s
#endif

    vendor_     = "(unknown)"
    version_    = "(unknown)"
    major_      = -1
    minor_      = -1
    patchlevel_ = -1

#if   defined (__GFORTRAN__)
    vendor_     = "GNU Fortran (GCC)"
    version_    = __VERSION__
    major_      = __GNUC__
    minor_      = __GNUC_MINOR__
    patchlevel_ = __GNUC_PATCHLEVEL__
#elif defined (__G95__)
    vendor_     = "G95 (GCC)"
    major_      = __G95__
    minor_      = __G95_MINOR__
    patchlevel_ = __G95_BUILD__
    write(version_,'(a,1x,i0,".",i0,1x,i0)') "g95", major_, minor_, patchlevel_
#elif defined (_CRAYFTN)
    vendor_     = "Cray Fortran (Crayftn)"
    version_    = _RELEASE_STRING
    major_      = _RELEASE_MAJOR
    minor_      = _RELEASE_MINOR
    i           = index (version_, "Version")
    if (i > 0) then
       !-----------------------------------------
       ! Parse version number from release string
       !-----------------------------------------
       s = version_(i+8:)
       forall(i=1:len(s)) s(i:i) = merge(" ",s(i:i),s(i:i)==".")
       read(s,*,iostat=ios) i, j, patchlevel_
       if (ios /= 0) patchlevel_ = -1
    end if
#elif defined (__INTEL_COMPILER)
    vendor_     = "Intel Fortran"
    i           = __INTEL_COMPILER
    select case (i)
    case (:1999)                    ! Older "Classic" Intel Compilers
       major_      =      i / 100
       i           = mod (i,  100)
       minor_      =      i / 10
#ifdef __INTEL_COMPILER_UPDATE  /* only for recent Intel version (>= 14?) */
       patchlevel_ = __INTEL_COMPILER_UPDATE
#else
       patchlevel_ = 0
#endif
    case (2000:9999)                ! OneAPI "Classic" Intel Compilers
       major_      = i
#if (__INTEL_COMPILER_UPDATE > 0) && (__INTEL_COMPILER_BUILD_DATE < 20240703)
       ! /* for several versions before deprecation except last one */
       minor_      = __INTEL_COMPILER_UPDATE
       patchlevel_ = 0
#else                           /* oneAPI deprecated versions... */
       version_    = compiler_version()
       i           = index (version_, "Version")
       if (i > 0) then
          !-----------------------------------------
          ! Parse version number from release string
          !-----------------------------------------
          s = version_(i+8:)
          forall(i=1:len(s)) s(i:i) = merge(" ",s(i:i),s(i:i)==".")
          read(s,*,iostat=ios) i, minor_, patchlevel_
          if (ios /= 0) minor_      = -1
          if (ios /= 0) patchlevel_ = -1
       end if
#endif
    case default                    ! OneAPI flang based Intel Compilers
       major_      =      i / 10000
       i           = mod (i,  10000)
       minor_      =      i / 100
       patchlevel_ = mod (i,  100)
    end select
    write(version_,'(i0,".",i0,".",i0,1x,"Build",1x,i0)') &
         major_, minor_, patchlevel_, __INTEL_COMPILER_BUILD_DATE
#elif defined (__ibm__)
    vendor_     = "IBM Fortran"
#elif defined (NAGFOR)
    vendor_     = "NAG Fortran"
#ifdef __NAG_COMPILER_RELEASE
    i           = __NAG_COMPILER_RELEASE    ! major*10+minor
    major_      =      i / 10
    minor_      = mod (i,  10)
    patchlevel_ = __NAG_COMPILER_BUILD
#endif
#elif defined (__NEC__)
    vendor_     = "NEC Fortran"
    version_    = compiler_version()
    !---------------------------------------------------
    ! Work around bad ASCII-NUL character in nfort-2.5.1
    !---------------------------------------------------
    s           = version_
    j           = min (len (version_), len (s))
    forall(i=1:j) version_(i:i) = merge(" ",s(i:i),s(i:i)==achar(0))
    i           = index (version_, ")")
    if (i > 0) then
       s = version_(i+2:)
       forall(i=1:len(s)) s(i:i) = merge(" ",s(i:i),s(i:i)==".")
       read(s,*,iostat=ios) i, j, patchlevel_
       if (ios == 0) then
          major_ = i
          minor_ = j
       else
          patchlevel_ = -1
       end if
    end if
#elif defined (__PGI)
    vendor_     = "PGI (pgfortran)"
#ifdef __PGIC__
    major_      = __PGIC__
    minor_      = __PGIC_MINOR__
#endif
#ifdef __PGIC_PATCHLEVEL__
    patchlevel_ = __PGIC_PATCHLEVEL__
    write(version_,'(i0,".",i0,".",i0)') major_, minor_, patchlevel_
#else
    write(version_,'(i0,".",i0)')        major_, minor_
#endif
#elif defined (__FLANG)
    vendor_     = "Flang (flang)"
    version_    = compiler_version ()
    i           = index (version_, "- ")
    if (i > 0) then
       !-------------------------------------------
       ! Parse version number from compiler_version
       !-------------------------------------------
       s = version_(i+2:)
       version_ = s
       forall (i=1:len(s)) s(i:i) = merge(" ",s(i:i),s(i:i)==".")
       read (s,*,iostat=ios) i, j
       major_ = i
       minor_ = j
    end if
#elif defined (__SUNPRO_F95)
    vendor_     = "Oracle Solarisstudio Fortran"
    write (s,95)
95  format (8h __SUNPRO_F95       )     ! Careful: Hex -> Hollerith!
    call hex2int (s, i)
    major_      =      i / 256
    minor_      = mod (i / 16,  16)
    patchlevel_ = mod (i,       16)
!   version_    = adjustl (s)
    write(version_,'(i0,".",i0,1x,"(",a,")")') &
         major_, minor_, trim (adjustl (s))
#endif
    if (present (vendor    )) vendor       = vendor_
    if (present (version   )) version      = version_
    if (present (major     )) major        = major_
    if (present (minor     )) minor        = minor_
    if (present (patchlevel)) patchlevel   = patchlevel_
  end subroutine compiler_info
  !---------------------------------------------------------------------
  subroutine hex2int (s, n)
    character(len=*), intent(in)  :: s
    integer,          intent(out) :: n
    integer :: i, k
    character :: c
    n = 0
    i = index (s, "0x")
    if (i > 0) then
       i = i + 2
       do while (i <= len_trim (s))
          c = s(i:i)
          read(c,'(z1)',err=999) k
          n = n * 16 + k
          i = i + 1
       end do
999    continue
    end if
  end subroutine hex2int
!------------------------------------------------------------------------------
end module mo_version
