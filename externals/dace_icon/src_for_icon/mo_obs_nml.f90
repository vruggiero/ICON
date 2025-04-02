!
!+ Read namelist /DEF_OBS_NML/
!
MODULE mo_obs_nml
!
! Description:
!   Read namelist /DEF_OBS_NML/.
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
!  remove old obsolete whitelist functionality
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup: remove unused variables
! V1_23        2013-03-26 Robin Faulwetter
!  Removed old reading routines for TOVS data (1dvar and sat_pp)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2008  original source
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp, sp           ! working precision kind parameter
  use mo_exception,  only: finish,         &! error exit routine
                           message          ! write warning message
  use mo_namelist,   only: position_nml,   &! position namelist
                           nnml,           &! namelist Fortran unit number
                           POSITIONED       ! ok    code from position_nml
  use mo_run_params, only: p_readbufr,     &! PE used to read BUFR data
                           ana_time         ! analysis time
  use mo_t_datum,    only: set_datum        ! routine to set datum
! use mo_t_use,      only: STAT_PASSIVE     !
  use mo_t_obs,      only: t_obs,          &! observation data type
                           t_spot,         &! observation meta data type
                           OBS_H,          &! observation id
                           OBS_T,          &!
                           OBS_U,          &!
                           OBS_V,          &!
                           TEMP             !
  use mo_temp,       only: t_temp,         &! TEMP level data type
                           check_store_temp ! check and store TEMP info
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_bcast          ! broadcast routine
  use mo_fdbk_tables,only: OT_TEMP          !
  use mo_obs_tables, only: check_report_0, &!
                           check_report_1
  implicit none
  !----------------
  ! public entities
  !----------------
  private
  public :: read_nml_obs

contains
!------------------------------------------------------------------------------
  subroutine read_nml_obs (obs)

  type (t_obs)  ,intent(inout) :: obs  ! observations data type to set
    !---------------------
    ! Namelist DEF_OBS_NML
    !---------------------
    real(wp)         :: lon    ! latitude  [degree]
    real(wp)         :: lat    ! longitude [degree]
    real(wp)         :: p      ! pressure  [hPa]
    character(len=8) :: rtype  ! report type
    character(len=2) :: otype  ! 'h', 't', 'rh', 'u', 'v'
    real(wp)         :: value  ! gpm   K    1    m/s  m/s
    integer          :: nlat   ! number of latitudes
!!! integer          :: prcflg ! processing flag
    namelist  /DEF_OBS_NML/ lon, lat, p, otype, rtype, value, nlat
    !----------------
    ! local variables
    !----------------
    logical                :: first
    integer                :: ierr
    integer                :: ident
    type (t_spot)          :: empty_spot
    type (t_spot)          :: spt
    type (t_temp)          :: empty_temp
    type (t_temp) ,pointer :: tmp(:)
    integer                :: novalid
    real(sp)               :: p_sp, val_sp
#if defined(__ibm__)
    integer                :: ios
#endif
    !---------------------------------
    ! loop over namelist group entries
    !---------------------------------
    first  = .true.
    ident  = 99000
    do
      !-------------
      ! set defaults
      !-------------
      lon    =   0._wp
      lat    =  50._wp
      p      = 500._wp ! hPa
      otype  = ''
      rtype  = 'TEMP'
      value  = 0._wp   ! K
      nlat   = 0
!!!   prcflg = 0
      !--------------
      ! read namelist
      !--------------
      if (dace% lpio) then
        call position_nml ('DEF_OBS_NML', lrewind=first, status=ierr)
        first = .false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=DEF_OBS_NML, iostat=ios)
          if (ios/=0) call finish ('read_nml_obs',                  &
                                   'ERROR in namelist /DEF_OBS_NML/')
#else
          read (nnml ,nml=DEF_OBS_NML)
#endif
        end select
      endif
      !------------------------------------------------------
      ! broadcast quantities, exit if namelist is not present
      !------------------------------------------------------
      call p_bcast (ierr   ,dace% pio)
      call p_bcast (lon    ,dace% pio)
      call p_bcast (lat    ,dace% pio)
      call p_bcast (p      ,dace% pio)
      call p_bcast (otype  ,dace% pio)
      call p_bcast (rtype  ,dace% pio)
      call p_bcast (value  ,dace% pio)
      call p_bcast (nlat   ,dace% pio)
!!!   call p_bcast (prcflg ,dace% pio)
      if (ierr/=POSITIONED) exit
      !FIXME
      call finish('read_nml_obs','FIXME: module currently broken')
      !------------------------------------------------
      ! set meta information, skip observation on error
      !------------------------------------------------
      if (dace% pe==p_readbufr) then
        select case (rtype)
        case default
          call message ('read_nml_obs',&
                        'namelist /DEF_OBS_NML/: invalid rtype = '//rtype)
        case ('TEMP')
          allocate (tmp (1))
          p_sp   = p
          val_sp = value
          ident  = ident + 1
          tmp               = empty_temp
          spt               = empty_spot
          spt% hd% modtype  = TEMP
          spt% hd% obstype  = OT_TEMP
          spt% hd% dbkz     = 520
!         spt% buf_type     =
!         spt% buf_subtype  =
          spt% col% nlev    = 1
          spt% col% c% dlat = lat
          spt% col% c% dlon = lon
!FIXME        spt% col% np      = 1
          spt% z            =      0._wp
          spt% ps           = 100000._wp
!         spt% stname       =
          spt% ident        = ident
          select case (otype)
          case ('h')
!FIXME            spt% col% ipar = OBS_H
            call set_datum (tmp(1)% gp, val_sp, 0)
!!!         tmp(1)% gp% prc = prcflg
          case ('t')
!FIXME            spt% col% ipar = OBS_T
            call set_datum (tmp(1)% t, val_sp, 0)
!!!         tmp(1)% t% prc = prcflg
          case ('u')
!FIXME            spt% col% ipar = OBS_U
            call set_datum (tmp(1)% uu, val_sp, 0)
!!!         tmp(1)% uu% prc = prcflg
          case ('v')
!FIXME            spt% col% ipar = OBS_V
            call set_datum (tmp(1)% vv, val_sp, 0)
!!!         tmp(1)% vv% prc = prcflg
!         case ('rh')
!           spt% col% ipar = OBS_RH
!           call set_datum (tmp(1)% rh, val_sp, 0)
!           tmp(1)% rh% prc = prcflg
          case default
            cycle
          end select
          call set_datum (tmp(1)% p, p_sp * 100._sp, 0)
          !---------
          ! set time
          !---------
          spt% hd% time        = ana_time
          spt%     actual_time = ana_time
!         call init_time (spt% actual_time, yyyy, mo, dd, hh, mi)
          !---------------------------
          ! insert in observation list
          !---------------------------
!         call check_report_0   (spt% use, spt% hd, 1, white=STAT_PASSIVE)
          call check_report_0   (spt% use, spt% hd, 1)
          call check_report_1   (spt)
          call check_store_temp (tmp, spt, obs, novalid)
          deallocate (tmp)
!         case ('RAD','TOVS')
!           art_nlon = nlat * 2
        end select
      endif
    end do
  end subroutine read_nml_obs

end module mo_obs_nml
