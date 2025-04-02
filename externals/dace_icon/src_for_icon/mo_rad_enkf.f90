!
!+ SEVIRI observation handling for LETKF (may be extended for other
!     instruments in future)
!
MODULE mo_rad_enkf
!
! Description:
!   SEVIRI observation handling
!
! Current Maintainer: DWD, Hendrik Reich, Lilo Bach, Annika Schomburg
!    phone: +49 69 8062 4943
!    fax:   +49 69 8062 3721
!    email: hendrik.reich@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_28        2014/02/26 Andreas Rhodin
!  new module for SEVIRI assimilation in LETKF
! V1_31        2014-08-21 Andreas Rhodin
!  disable temporary fixes for SEVIRI
! V1_45        2015-12-15 Andreas Rhodin
!  pass H_x0 (ensemble mean) to check_seviri
! V1_51        2017-02-24 Axel Hutt
!  namelist extension in mo_seviri
!              2022-09-22 Annika Schomburg
!  clean up
!              2022-09-20 Annika Schomburg
! rename: from mo_seviri.f90 to mo_rad_enkf.f90
!              2022-09-30 Annika Schomburg
! new subroutine: seviri_wv_error_model
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

  !=============
  ! modules used
  !=============
  use mo_exception,   only: finish         ! abort routine
  use mo_mpi_dace,    only: dace,         &! MPI group info
                            p_bcast        ! broadcast routine
  use mo_namelist,    only: position_nml, &! position namelist
                            nnml,         &! namelist Fortran unit number
                            POSITIONED     ! ok    code from position_nml

  use mo_kind,         only: wp            ! working precision kind
  use mo_instrid,      only: rttov_instr, &! RTTOV from WMO instrument id
                              hss_instr

  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: read_nml_seviri   ! read namelist /SEVIRI_OBS/
  public :: read_letkf_obserr_rad_nml! read namelist /LETKF_OBSERR_RAD/
  public :: seviri_wv_error_model ! subroutine computing obs errors
  public :: hyp_ir_error_model! computes hyperspectral observation error
  public :: seviri_error_model! Observation error model
  public :: plevel_select     ! determine plevel from ensemble:
                              ! 0=ensemble mean, 1=choose highest plevel
  public :: letkf_obserr_ir   ! stores obserr info for different IR instruments
  public :: n_letkf_obserr_rad! number of entries in letkf_obserr_ir
!------------------------------------------------------------------------------
  !=================
  ! module variables
  !=================
  type t_letkf_obserr_rad
    integer  :: verbose            = 1
    integer  :: instr              = -1
    integer  :: rad_error_model    = -1
    real(wp) :: ens_cld_limit      = -1
  end type t_letkf_obserr_rad


  integer     :: n_letkf_obserr_rad = 0

  type(t_letkf_obserr_rad), target :: letkf_obserr_ir(5)


  !----------------------
  ! Namelist /SEVIRI_OBS/
  !----------------------
  integer  :: verbose            = 0        ! verbosity flag
  integer  :: seviri_error_model = 0        ! default: constant observation error
  integer  :: plevel_select      = 0        ! default: average plevels

  namelist /SEVIRI_OBS/ verbose, seviri_error_model, plevel_select

!==============================================================================
contains
!==============================================================================

  subroutine read_nml_seviri
  !---------------------------
  ! read namelist /SEVIRI_OBS/
  !---------------------------
    integer :: ierr
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('SEVIRI_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=SEVIRI_OBS, iostat=ierr)
        if(ierr/=0) call finish('read_nml_seviri','ERROR in namelist /SEVIRI_OBS/')
#else
        read (nnml ,nml=SEVIRI_OBS)
#endif
      end select
    endif
    !------------------------
    ! broadcast to other PE's
    !------------------------
    call p_bcast (verbose           ,dace% pio)
    call p_bcast (seviri_error_model,dace% pio)
    call p_bcast  (plevel_select,     dace% pio)
    !----------------------------
    ! print namelist /SEVIRI_OBS/
    !----------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  namelist /SEVIRI_OBS/'
      write(6,'()')
      write(6,'(a,i6,4x,a)') 'verbose            =', verbose ,'verbosity flag'
      write(6,'(a,i6     )') 'seviri_error_model =', seviri_error_model
      write(6,'(a,i6     )') 'plevel_select      =', plevel_select
      write(6,'()')
    endif
    !---------------------------------
    ! check namelist for plausibility
    !---------------------------------
    if ( seviri_error_model < 0 .or. seviri_error_model > 2 ) then
       call finish("read_nml_seviri", "invalid seviri_error_model")
    end if
    if ( plevel_select < 0 .or. plevel_select > 1 ) then
       call finish("read_nml_seviri", "invalid plevel_select")
    end if
  end subroutine read_nml_seviri



!--------------------------------------------------------------------------------
  subroutine seviri_wv_error_model(channel,cl_impact,eo)
! This subroutine computes the observation errors for the allsky assimilation of
! the SEVIRI WV channels based on a polynomial fit of the standard deviation of
! the first guess departures with the cloud impact (=symmetric difference
! between allsky and clear-sky brightness temperatures) as predictor
!
   integer, intent(in)            :: channel       ! satellite channel
   real(wp), intent(in)           :: cl_impact   ! cloud_impact
   real(wp), intent(out)          :: eo       ! observation error

   real(wp)               :: tbdiff_max, err_max, bcoeff(4), px, infl
   real(wp)               :: minerr, fitconst


   minerr = 2.0_wp  !error should not be lower than this value
   infl = 1.5_wp
   if (channel == 5) then   !SEVIRI channel 5
     tbdiff_max = 9.0_wp
     err_max = 5.3258_wp
     bcoeff = (/ 3.0953, 0.7701, -0.0116, -0.0103 /)
     fitconst = 4.5_wp
     px = cl_impact -fitconst
   else if (channel == 6)  then   !SEVIRI channel 6
     tbdiff_max = 18.0_wp
     err_max = 9.8141
     bcoeff = (/5.8082, 0.619, -0.0086, -0.0009 /)
     fitconst = 9.0_wp
     px = cl_impact -fitconst
   end if
   if (cl_impact  >= tbdiff_max ) then
     eo = err_max + infl
   else
     !eo = bcoeff(1) + bcoeff(2) * px + bcoeff(3)*px**2 + bcoeff(4)*px**3 + infl
     !Horner-Scheme applied to equation above (more efficient):
     eo = bcoeff(1) + (bcoeff(2) + (bcoeff(3) + bcoeff(4)*px) *px)*px + infl
   end if
   eo = MAX(eo,minerr)

  end subroutine seviri_wv_error_model
!----------------------------------------------------------------------------------

  subroutine read_letkf_obserr_rad_nml
  !=================================
  ! read namelist /LETKF_OBSERR_RAD/
  !=================================

    !----------------------------
    ! Namelist /LETKF_OBSERR_RAD/
    !----------------------------
    integer  :: instr                         ! WMO instrument Id
    integer  :: rad_error_model               ! default: constant observation error
    integer  :: verbose                       ! verbosity flag
    real(wp) :: ens_cld_limit                 ! maximum cloudiness-agreement allowed among ensemble members

    namelist /LETKF_OBSERR_RAD/ instr, rad_error_model, verbose, ens_cld_limit
    !----------------
    ! local variables
    !----------------
    type(t_letkf_obserr_rad), pointer :: oep => null()
    logical                           :: first
    integer                           :: ierr
    integer                           :: sensor ! RTTOV instrument ID
#if defined(__ibm__)
    integer                           :: ios
#endif

    first = .true.
    do
      !-------------
      ! set defaults
      !-------------
      instr                = -1
      rad_error_model      = 0
      verbose              = 0
      ens_cld_limit        = 0.2_wp
      !----------------------------------
      ! read namelist, consistency check
      !----------------------------------
      if (dace% lpio) then
        call position_nml('LETKF_OBSERR_RAD', lrewind=first, status=ierr)
        first=.false.
        select case (ierr)
        case(POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=LETKF_OBSERR_RAD, iostat=ios)
          if(ios/=0) call finish('read_letkf_obserr_rad_nml','ERROR in namelist /LETKF_OBSERR_RAD/')
#else
          read (nnml ,nml=LETKF_OBSERR_RAD)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (instr,             dace% pio)
      call p_bcast (rad_error_model,   dace% pio)
      call p_bcast (verbose,           dace% pio)
      call p_bcast (ens_cld_limit,     dace% pio)
      !----------------
      n_letkf_obserr_rad = n_letkf_obserr_rad + 1
      if (n_letkf_obserr_rad > size(letkf_obserr_ir)) &
      call finish('read_letkf_obserr_rad_nml', 'n_letkf_obserr_rad > 5')

      oep => letkf_obserr_ir(n_letkf_obserr_rad)

      oep% instr            = instr
      oep% rad_error_model  = rad_error_model
      oep% verbose          = verbose
      oep% ens_cld_limit    = ens_cld_limit
      !----------------------------
      ! print namelist /LETKF_OBSERR_RAD/
      !----------------------------
      if (dace% lpio) then
        write(6,'(a)') repeat('-',79)
        write(6,'()')
        write(6,'(a)') '  namelist /LETKF_OBSERR_RAD/'
        write(6,'()')
        write(6,'(a,i6     )') 'instr              =', instr
        write(6,'(a,i6     )') 'rad_error_model    =', rad_error_model
        write(6,'(a,f6.2  )')  'ens_cld_limit      =', ens_cld_limit
        write(6,'(a,i6,4x,a)') 'verbose            =', verbose ,'verbosity flag'
        write(6,'()')
      endif
      !---------------------------------
      ! check namelist for plausibility
      !---------------------------------
      if ( rad_error_model < 0 .or. rad_error_model > 2 ) then
        call finish("read_letkf_obserr_rad_nml", "invalid rad_error_model")
      end if

      sensor = rttov_instr(instr) ! RTTOV ID
      if (.not. hss_instr(sensor)) then ! if not hyperspectral sounders
        call finish("read_letkf_obserr_rad_nml", "invalid WMO instr ID")
      end if

      if ( ens_cld_limit < 0 .or. ens_cld_limit > 1 ) then
        call finish("read_letkf_obserr_rad_nml", "invalid ens_cld_limit")
      end if
    end do

  end subroutine read_letkf_obserr_rad_nml

subroutine hyp_ir_error_model(channel, x, eo, limit)
  integer,  intent(in)    :: channel
  real(wp), intent(in)    :: x !ensmean_cld_fl
  real(wp), intent(in)    :: limit
  real(wp), intent(inout) :: eo

  real(wp)                :: minerr, infl

  minerr = eo
  infl   = 1.5_wp

  if (x > limit ) then
    eo = 90.0_wp
  else
    eo = minerr + (x * infl)
  end if
end subroutine hyp_ir_error_model


!==============================================================================


end module mo_rad_enkf
