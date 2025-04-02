!
!+ Write profiles of atmospheric parameters (for gnuplot)
!
MODULE mo_profile
!
! Description:
!   Routines to write profiles of atmospheric parameters
!   (suitable for plotting with gnuplot).
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
! V1_7         2009/08/24 Harald Anlauf
!  write_profile : use correct number of levels for isobaric grid
! V1_8         2009/12/09 Andreas Rhodin
!  work around bug in xlf V12.1 (internal compiler error)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Reduce debugging output
! V1_20        2012-06-18 Andreas Rhodin
!  rename mo_profiles to mo_profile (make module name and file name consistent)
! V1_37        2014-12-23 Harald Anlauf
!  print additional information for nearest neighbors; option vor area averaging
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Autors:
! A.Rhodin  DWD  2007  original source
!------------------------------------------------------------
!-------------
! modules used
!-------------
#if defined(__ibm__)
use mo_exception,    only: finish            ! exit on error condition
#endif
use mo_kind,         only: wp                ! working precision kind parameter
use mo_mpi_dace,     only: dace,            &! MPI group info
                           p_bcast,         &! generic MPI broadcast routine
                           p_sum             ! sum over PEs
use mo_namelist,     only: position_nml,    &! routine to position nml group
                           nnml,            &! namelist fortran unit number
                           POSITIONED        ! position_nml: OK return flag
use mo_atm_state,    only: t_atm             ! atmospheric state derived type
use mo_wmo_tables,   only: WMO3_ISOBARIC,   &! level type
                           WMO6_LATLON,     &! grid types
                           WMO6_ROTLL,      &!
                           WMO6_GAUSSIAN     !
!use coordinates,     only: cartesian         !+++ required for IBM clf V12.1
use mo_grid_intpol,  only: Grid_Indices      ! get indices of neighbor gridpts.
use mo_fortran_units,only: get_unit_number,& ! get unused unit number
                           return_unit_number! release unit number after use
use mo_run_params   ,only: aux,             &! path for auxiliary output
                           path_file         ! concatenate path+filename
use mo_physics      ,only: d2r,             &! factor: degree -> radians
                           r2d,             &! factor: radians -> degree
                           pi                ! 3.1415...
implicit none

!----------------
! public entities
!----------------
private
public :: read_nml_profiles  ! read namelist /PROFILES/
public :: write_profiles     ! write profiles at points specified by namelist
public :: write_profile      ! write an atmospheric profile

!---------
! namelist
!---------
integer ,parameter :: mp = 20         ! max no. of profiles

real(wp) :: lat_lon (2,mp) = -999._wp ! coordinates of profiles
integer  :: ih             =    2     ! 1=next, 2=interpolation, 3=area average
namelist /PROFILES/ lat_lon, ih

contains
!------------------------------------------------------------------------------
  subroutine read_nml_profiles
  !-------------------------
  ! read namelist /PROFILES/
  !-------------------------
    integer :: ierr, i
#if defined(__ibm__)
    integer :: ios
#endif
    lat_lon = -999._wp ! coordinates of profiles
    ih      =    2     ! 1=next, 2=interpolation, 3=area average
    if (dace% lpio) then
      call position_nml ('PROFILES' ,lrewind=.true. ,status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=PROFILES, iostat=ios)
        if (ios/=0) call finish ('read_nml_profiles',&
                                 'ERROR in namelist /PROFILES/')
#else
        read (nnml ,nml=PROFILES)
#endif
      end select
    endif
    call p_bcast (lat_lon, dace% pio)
    call p_bcast (ih,      dace% pio)
    !--------------------------
    ! print namelist /PROFILES/
    !--------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  namelist /PROFILES/'
      write(6,'()')
      select case (ih)
      case (1)
        write(6,'(a,i16,a)') '    ih      = ',ih,' (nearest neighbour)'
      case (2)
        write(6,'(a,i16,a)') '    ih      = ',ih,' (interpolation)'
      case (3)
        write(6,'(a,i16,a)') '    ih      = ',ih,' (area averaging)'
      end select
      do i = 1, mp
        if (abs(lat_lon(1,i)) >  90._wp) exit
        if (abs(lat_lon(2,i)) > 360._wp) exit
        write(6,'(a,2f8.3)') '    lat,lon = ',lat_lon(:,i)
        if (ih == 3 .and. mod(i,2)==0) write(6,'(a)')
      end do
      write(6,'(a)')
    endif

  end subroutine read_nml_profiles
!------------------------------------------------------------------------------
  subroutine write_profiles (atm, name, comment)
  type (t_atm)     ,intent(in) :: atm     ! atmospheric state to plot
  character(len=*) ,intent(in) :: name    ! prefix of file name
  character(len=*) ,intent(in) :: comment ! comment line to write to file
                    optional   :: comment
  !---------------------------------------------------
  ! write profiles at points specified by the namelist
  !---------------------------------------------------
    integer          :: ip    ! profile index
    character(len=3) :: suff  ! suffix of file name
    !-----------------
    ! write to logfile
    !-----------------
    if (dace% lpio) then
!     write (6,'(a)') repeat('-',79)
      write (6,'( )')
      write (6,'(a,a)') '  writing atmospheric profiles: ',trim(name)
      write (6,'( )')
    endif

    select case (ih)
    case (3)
      !-------------------------
      ! loop over pairs of edges
      !-------------------------
      do ip = 2, mp, 2
        write (suff,'("_",i2.2)') ip / 2
        call write_profile (atm, trim(name)//suff, ih,        &
                            lat_lon(1,ip-1), lat_lon(2,ip-1), &
                            lat_lon(1,ip  ), lat_lon(2,ip  ), &
                            comment=comment)
      end do
    case default
      !-------------------------
      ! loop over profile points
      !-------------------------
      do ip = 1, mp
        write (suff,'("_",i2.2)') ip
        call write_profile (atm, trim(name)//suff, ih,     &
                            lat_lon(1,ip), lat_lon(2,ip),  &
                            comment=comment)
      end do
    end select
    if (dace% lpio) write (6,'( )')
  end subroutine write_profiles
!------------------------------------------------------------------------------
  subroutine write_profile (atm, name, ih, lat, lon, lat2, lon2 , comment)
  type (t_atm)     ,intent(in) :: atm     ! atmospheric state to plot
  character(len=*) ,intent(in) :: name    ! file basename
  integer          ,intent(in) :: ih      ! interpolation flag
  real(wp)         ,intent(in) :: lat     ! latitude
  real(wp)         ,intent(in) :: lon     ! longitude
  real(wp)         ,intent(in) :: lat2    ! latitude  of 2nd edge (if ih==3)
  real(wp)         ,intent(in) :: lon2    ! longitude of 2nd edge (if ih==3)
  character(len=*) ,intent(in) :: comment ! comment line to write to file
                    optional   :: comment, lat2, lon2
  !-----------------------------
  ! write an atmospheric profile
  !-----------------------------
    integer    ,parameter :: map = 20     ! max.no. atmospheric parameters
    integer    ,parameter :: ngp = 16     ! max.no. interpolation points
    real(wp)              :: w   (ngp)    ! weight
    real(wp) ,allocatable :: wg (:,:,:)   ! grid point weights
    integer               :: idx (ngp,4)  ! indices
    integer               :: np           ! no.interpolation points used
    integer               :: i            ! index (atmospheric parameter)
    integer               :: j            ! index (interpolation coefficient)
    integer               :: k            ! index (level)
    integer               :: i1,i2,i4     ! model grid indices
    real(wp)              :: rlon, rlat   ! edges (radiand)
    real(wp)              :: rlon2, rlat2 ! edges (radiand)
    real(wp)              :: glon         ! longitude of gridpoint
    real(wp) ,allocatable :: tmp (:,:)    ! temporary
    integer               :: nz           ! number of layers
    integer               :: nap          ! number of atmospheric parameters
    character(len=16)     :: names(map)   ! names of atmospheric parameters
    integer               :: iu           ! Fortran unit number
    character(len=128)    :: filename
    integer               :: nnx (4)      ! nearest neighbor indices

    nnx = 0
    !----------------------------
    ! check for valid coordinates
    !----------------------------
    if (abs(lat)>  90._wp) return
    if (abs(lon)> 360._wp) return
    if (dace% lpio) write (6,'(a,2f10.3)') '    lat, lon   : ',lat,lon
    select case (ih)
    case (3)
      !------------------
      ! determine weights
      !------------------
      if (dace% lpio) write (6,'(a,2f10.3)') '    lat, lon   : ',lat2,lon2
      rlat  = d2r * lat
      rlon  = d2r * lon
      rlat2 = d2r * lat2
      rlon2 = d2r * lon2
      allocate (wg (atm% lb (1): atm% ub (1) &
                   ,atm% lb (2): atm% ub (2) &
                   ,atm% lb (4): atm% ub (4)))
      wg   = 0._wp
      w(1) = 1._wp
      do i4 = atm% lb (4), atm% ub (4)
        do i2 = atm% lb (2), atm% ub (2)
          select case (atm% grid% gridtype)
          case (WMO6_LATLON, WMO6_ROTLL, WMO6_GAUSSIAN)
            w(1) = cos (d2r * atm% grid% dlat (i2))
          end select
          do i1 = atm% lb (1), atm% ub (1)
            if (atm% grid% marr (1,i1,i2,i4) /= dace% pe) cycle
            if (atm% grid% rlat (i1,i2,1,i4) <  rlat    ) cycle
            if (atm% grid% rlat (i1,i2,1,i4) >  rlat2   ) cycle
            glon = atm% grid% rlon (i1,i2,1,i4)
            if (glon < rlon ) glon = glon + 2._wp * pi
            if (glon > rlon2) glon = glon - 2._wp * pi
            if (glon < rlon ) cycle
            wg (i1,i2,i4) = w(1)
          end do
        end do
      end do
      w(1) = p_sum (sum (wg))
      if (w(1)==0.) then
        if (dace% lpio) &
          write(6,'(a,4f8.2)') '  no gridpoints in area ',lat,lat2,lon,lon2
        return
      endif
      wg = wg / w(1)
    case default
      !-------------------------------------
      ! determine interpolation coefficients
      !-------------------------------------
      call grid_indices &
        (lon,       & ! <-- geodetic longitude
         lat,       & ! <-- geodetic latitude
         atm% grid, & ! <-- grid data type
         idx,       & ! --> Grid point indices [Point, index]
         w,         & ! --> Weights
         np)          ! --> number of points returned
      !------------------------------------
      ! skip if coordinates are out of area
      !------------------------------------
      if (np==0) then
        if (dace% lpio) &
          write(6,'(a,2f8.2)') '  profile location out of area:',lat, lon
        return
      endif
      !-----------------------------------------------
      ! modify weights for nearest neighbour selection
      !-----------------------------------------------
      if (ih==1) then
        j    = sum(maxloc(w(1:np)))
        w    = 0._wp
        w(j) = 1._wp
        nnx  = idx(j,:)
        if (dace% lpio) then
          write(6,'(a,2f10.3)') '    nearest gridpoint:',      &
               atm% grid% rlat (nnx(1),nnx(2),1,nnx(3)) * r2d, &
               atm% grid% rlon (nnx(1),nnx(2),1,nnx(3)) * r2d
          write(6,'(a,3i10)')   '    with grid indices:', nnx(1),nnx(2),nnx(3)
        endif
      endif
    end select
    !---------------------
    ! allocate temporaries
    !---------------------
    select case (atm% grid% levtyp)
    case (WMO3_ISOBARIC)
       nz = atm% grid% nz
    case default
       nz = atm% grid% nz + 1
    end select
    allocate (tmp (nz, map))
    tmp (:,:) = 0._wp
    !---------------------------------
    ! loop over atmospheric parameters
    !---------------------------------
    nap = 0
    do i = 1, size(atm% m)
      !------------------------------
      ! pick up values to interpolate
      !------------------------------
      if (      atm% m(i)% i% ub(3) < nz-1 ) cycle
      if (.not. atm% m(i)% i% alloc        ) cycle
      if (      atm% m(i)% i% rep   /= 'gg') cycle
      nap = nap + 1
      names(nap) = atm% m(i)% i% name
      if (dace% lpio) then
        write (6,'(a,a)') '      parameter: ',names(nap)
      endif
      !-----------------
      ! loop over levels
      !-----------------
      tmp (:, nap) = sqrt(huge(1._wp))
      do k = atm% m(i)% i% lb(3), atm% m(i)% i% ub(3)
        if (k<1 .or. k>nz)  cycle
        select case (ih)
        case (3)
          !--------
          ! average
          !--------
          tmp (k, nap) = sum (wg * atm% m(i)% ptr (:,:,k,:))
        case default
          !------------
          ! interpolate
          !------------
          tmp (k, nap) = 0._wp
          do j = 1, np
            if (idx(j,4) /= dace% pe) cycle
            tmp (k, nap) = tmp (k, nap)                                       &
                         + w(j) * atm% m(i)% ptr (idx(j,1),idx(j,2),k,idx(j,3))
          end do
        end select
      end do
      if (nap==map) exit
    end do
    !---------------
    ! send to I/O PE
    !---------------
    tmp = p_sum (tmp)
    !---------------
    ! write plotfile
    !---------------
    if (dace% lpio) then
      iu = get_unit_number()
      filename = path_file(aux,trim(name)//'.prf')
      open  (iu, file=filename, action='write')
      write (iu,'(a           )')'#'
      write (iu,'(a           )')'# atmospheric profiles'
      write (iu,'(a,a         )')'# ',trim(name)
      if(present(comment)) &
        write (iu,'(a,a       )')'# ',trim(comment)
      write (iu,'(a           )')'#'
      select case (ih)
      case (1)
        write (iu,'(a,f10.3   )')'# latitude      =',lat
        write (iu,'(a,f10.3   )')'# longitude     =',lon
        write (iu,'(a,i6,a    )')'# interpolation =',ih,' (nearest neighbour)'
        write (iu,'(a,2f10.3  )')'#  nearest gridpoint:',    &
             atm% grid% rlat (nnx(1),nnx(2),1,nnx(3)) * r2d, &
             atm% grid% rlon (nnx(1),nnx(2),1,nnx(3)) * r2d
        write (iu,'(a,3i10    )')'#  with grid indices:', nnx(1),nnx(2),nnx(3)
      case (2)
        write (iu,'(a,f10.3   )')'# latitude      =',lat
        write (iu,'(a,f10.3   )')'# longitude     =',lon
        write (iu,'(a,i6,a    )')'# interpolation =',ih,' (interpolated)'
      case (3)
        write (iu,'(a,2f10.3  )')'# latitude      =',lat, lat2
        write (iu,'(a,2f10.3  )')'# longitude     =',lon, lon2
        write (iu,'(a,i6,a    )')'# interpolation =',ih,' (area average)'
      end select
      write (iu,'(a           )')'#'
      write (iu,'(a,i5,a,a    )')'# column',1,  '   :  index'
      do i = 1, nap
        write (iu,'(a,i5,a,a  )')'# column',i+1,'   : ',names(i)
      end do
      write (iu,'(a           )')'#'
      write (iu,'(a,20(1x,a16))')'# i        ',names(1:nap)
      write (iu,'(a           )')'#'
      do k=1,nz
        write (iu,'(i3,20(1x,f16.6))') k, (tmp (k,i), i=1,nap)
      end do
      close (iu)
      call return_unit_number(iu)
    endif
    !---------
    ! clean up
    !---------
    deallocate (tmp)
  end subroutine write_profile
!------------------------------------------------------------------------------
end module mo_profile
