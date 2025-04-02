!
!+ make available CADS v3.1 via a module
!
module mo_cads31
!
! Description:
! This module contains the very slightly modified CADS V3.1 code
! supplemented by routines, that are required for using CADS in
! the DWD code.
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2020/12/15 Robin Faulwetter
!  intial version
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Robin Faulwetter  DWD
! Olaf Stiller      DWD
!=======================================================================

  !=============
  ! Modules used
  !=============
  use mo_mpi_dace,          only: dace,  &
                                  p_bcast,&
                                  p_bcast_ptr
  use mo_exception,         only: finish
  use mo_kind,              only: wp, i4
  use utilities,            only: sortrx
  use CADS31_Module


  implicit none
  private

  public :: cads_setup_cloud_31
  public :: cads_detect_cloud_31
  public :: cads_setup_aerosol_31
  public :: cads_detect_aerosol_31
  public :: cads_setup_land_sensitivity_31
  public :: cads_detect_land_sensitivity_31
  public :: cads_setup_trace_gas_31
  public :: cads_detect_trace_gas_31
  public :: p_bcast_cloud_31
  public :: p_bcast_aerosol_31
  public :: p_bcast_land_sens_31
  public :: p_bcast_trace_gas_31
  public :: set_cloud_detect_chans_31
  public :: check_sensor_id_31
!  public :: get_cloud_opt_31
  public :: l_debug
  public :: dwd_sort
  public :: l_rank_p
  public :: t_aer_out
  public :: t_cld_out
  public :: usd

  type t_aer_out
    logical  :: l_out      = .false.
    integer  :: aer_type   = 0
    real(wp) :: hgt_thresh = 0._wp
    real(wp) :: aod_other  = 0._wp
    real(wp) :: aod_dust   = 0._wp
    real(wp) :: aod_ash    = 0._wp
  end type t_aer_out

  integer, parameter :: mx_cl = 7
  type t_cld_out
    logical  :: l_out                 = .false.
    integer  :: i__chan_low           = 0
    integer  :: i__chan_high          = 0
    integer  :: i__scenario_index     = 0
    integer  :: i__start_channel      = 0
    real(wp) :: z__cloud_level        = 0._wp
    real(wp) :: z__Wsqdev             = 0._wp
    real(wp) :: z__sqdev(mx_cl)       = 0._wp
    real(wp) :: z__intercluster_max   = 0._wp
    real(wp) :: z__intercluster_sqdev = 0._wp
  end type t_cld_out
 
    
  
  logical :: l_debug  = .false.
  logical :: dwd_sort = .false.
  logical :: l_rank_p = .false.
  integer :: usd      = 0

  interface cads_setup_cloud_31             ; module procedure cads_setup_cloud             ; end interface
  interface cads_detect_cloud_31            ; module procedure cads_detect_cloud            ; end interface
  interface cads_setup_aerosol_31           ; module procedure cads_setup_aerosol           ; end interface
  interface cads_detect_aerosol_31          ; module procedure cads_detect_aerosol          ; end interface
  interface cads_setup_land_sensitivity_31  ; module procedure cads_setup_land_sensitivity  ; end interface
  interface cads_detect_land_sensitivity_31 ; module procedure cads_detect_land_sensitivity ; end interface
  interface cads_setup_trace_gas_31         ; module procedure cads_setup_trace_gas         ; end interface
  interface cads_detect_trace_gas_31        ; module procedure cads_detect_trace_gas        ; end interface
  interface p_bcast_cloud_31                ; module procedure p_bcast_cloud                ; end interface
  interface p_bcast_aerosol_31              ; module procedure p_bcast_aerosol              ; end interface
  interface p_bcast_land_sens_31            ; module procedure p_bcast_land_sens            ; end interface
  interface p_bcast_trace_gas_31            ; module procedure p_bcast_trace_gas            ; end interface
  interface set_cloud_detect_chans_31       ; module procedure set_cloud_detect_chans       ; end interface
  interface check_sensor_id_31              ; module procedure check_sensor_id              ; end interface
!  interface get_cloud_opt_31                ; module procedure get_cloud_opt                ; end interface

contains

  ! Basics
  include "find_nml.incf"
  include "CADS31_Abort.incf"
  ! Cloud detection
  include "CADS31_Setup_Cloud.incf"
  include "CADS31_Detect_Cloud.incf"
  include "CADS31_Detect_Cloud_Heapsort.incf"
  include "CADS31_Detect_Cloud_Imager.incf"
  include "CADS31_Detect_Cloud_Scenario.incf"
  include "CADS31_Detect_Cloud_Separator.incf"
  include "CADS31_Detect_Cloud_Smooth.incf"
  ! Aerosol detection
  include "CADS31_Setup_Aerosol.incf"
  include "CADS31_Detect_Aerosol.incf"
  ! Trace_Gas detection
  include "CADS31_Setup_Trace_Gas.incf"
  include "CADS31_Detect_Trace_Gas.incf"
  ! Land_Sensitivity detection
  include "CADS31_Setup_Land_Sensitivity.incf"
  include "CADS31_Detect_Land_Sensitivity.incf"
  ! normalized height correction by DWD:
  include "nhgt2p.incf"

  subroutine cads_sort(n, rdata, iindex)
    integer,       intent(in)    :: n
    real(kind=wp), intent(in)    :: rdata(:)
    integer,       intent(inout) :: iindex(:)
    if (dwd_sort) then
      call sortrx(n, rdata, iindex)
    else
      call CADS_Detect_Cloud_Heapsort(n, rdata, iindex)
    end if
  end subroutine cads_sort

  function check_sensor_id(id, typ) result(l)
    logical                                :: l
    integer,          intent(in)           :: id
    character(len=*), intent(in), optional :: typ

    l = (id >= JP__Min_Sensor_Index) .and. (id <= JP__Max_Sensor_Index)
    if (l .and. present(typ)) then
      select case(typ)
      case('cloud')
        l = any(s__cads_setup_cloud    (:)%m__sensor == id)
      case('aerosol')
        l = any(s__cads_setup_aerosol  (:)%m__sensor == id)
      case('land')
        l = any(s__cads_setup_land     (:)%m__sensor == id)
      case('trace_gas')
        l = any(s__cads_setup_trace_gas(:)%m__sensor == id)
      case default
        call finish('check_sensor_id','typ='//trim(typ)//' not implemented.')
      end select
    end if

  end function check_sensor_id


  subroutine set_cloud_detect_chans(isens, nbands, nchans, chans)
    integer, intent(in) :: isens
    integer, intent(in) :: nbands
    integer, intent(in) :: nchans(:)
    integer, intent(in) :: chans(:,:)

    integer :: mx_chans

    if ( isens < lbound(s__cads_setup_cloud,1) .or. &
         isens > ubound(s__cads_setup_cloud,1))     &
         call finish('set_cloud_detect_chans', 'Invalid sensor in set_cloud_detect_chans')

    if ( .not.associated(s__cads_setup_cloud(isens)%n__band_size) .or. &
         .not.associated(s__cads_setup_cloud(isens)%n__bands))         &
      call finish('set_cloud_detect_chans',       &
           'Attempt to modify settings for sensor, that was not initialized')

    mx_chans = maxval(nchans)
    if (s__cads_setup_cloud(isens)%n__num_bands < nbands) then
      deallocate(s__cads_setup_cloud(isens)%n__band_size)
      allocate(s__cads_setup_cloud(isens)%n__band_size(nbands))
    end if
    if ( s__cads_setup_cloud(isens)%n__num_bands < nbands .or. &
         size(s__cads_setup_cloud(isens)%n__bands,1) < mx_chans) then
      deallocate(s__cads_setup_cloud(isens)%n__bands)
      allocate(s__cads_setup_cloud(isens)%n__bands(mx_chans,nbands))
    end if

    s__cads_setup_cloud(isens)%n__num_bands = nbands
    s__cads_setup_cloud(isens)%n__band_size(1:nbands) = nchans(1:nbands)
    s__cads_setup_cloud(isens)%n__bands(1:mx_chans,1:nbands) = chans(1:mx_chans,1:nbands)

  end subroutine set_cloud_detect_chans

  ! subroutine get_cloud_opt(isens, n__bandtouse)
  !   integer,          intent(in)            :: isens
  !   integer, pointer, intent(out), optional :: n__bandtouse(:)

  !   if (present(n__bandtouse)) then
  !     if (associated(s__cads_setup_cloud(isens)%n__bandtouse)) then
  !       n__bandtouse => s__cads_setup_cloud(isens)%n__bandtouse
  !     else
  !       n__bandtouse => null()
  !     end if
  !   end if

  ! end subroutine get_cloud_opt

  !=============
  ! MPI routines
  !=============

#define BCAST_P(x) if (.not.dace%lpio) x=>null(); call p_bcast_ptr(x,pio)

  subroutine p_bcast_cloud(pio)
    integer,      intent(in)    :: pio
    type(cloud_detect_type), pointer :: s => null()
    integer :: i

    do i = lbound(s__cads_setup_cloud,1), ubound(s__cads_setup_cloud,1)
      s => s__cads_setup_cloud(i)
      call p_bcast_cloud_type(s, pio)
      BCAST_P(s%n__gradchkinterval      )
      BCAST_P(s%n__band_size            )
      BCAST_P(s%n__bands                )
      BCAST_P(s%n__window_width         )
      BCAST_P(s%n__window_bounds        )
      BCAST_P(s%n__bandtouse            )
      BCAST_P(s%r__bt_threshold         )
      BCAST_P(s%r__grad_threshold       )
      BCAST_P(s%r__window_grad_threshold)
      BCAST_P(s%n__imager_chans         )
      BCAST_P(s%r__stddev_threshold     )
    end do

  end subroutine p_bcast_cloud

  subroutine p_bcast_aerosol(pio)
    integer,      intent(in)    :: pio
    type(aerosol_detect_type), pointer :: s => null()
    integer :: i

    do i = lbound(s__cads_setup_aerosol,1), ubound(s__cads_setup_aerosol,1)
      s => s__cads_setup_aerosol(i)
      call p_bcast_aerosol_type(s, pio)
      BCAST_P(s%n__num_regression   )
      BCAST_P(s%n__num_aerosol_chans)
      BCAST_P(s%n__aerosol_chans    )
      BCAST_P(s%r__aerosol_tbd      )
      BCAST_P(s%r__coef_aod         )
    end do

  end subroutine p_bcast_aerosol

  subroutine p_bcast_land_sens(pio)
    integer,      intent(in)    :: pio
    integer :: i

    do i = lbound(s__cads_setup_land,1), ubound(s__cads_setup_land,1)
       call p_bcast_land_sens_type(s__cads_setup_land(i), pio)
    end do

  end subroutine p_bcast_land_sens

  subroutine p_bcast_trace_gas(pio)
    integer,      intent(in)    :: pio
    type(trace_gas_detect_type), pointer :: s => null()
    integer :: i

    do i = lbound(s__cads_setup_trace_gas,1), ubound(s__cads_setup_trace_gas,1)
      s => s__cads_setup_trace_gas(i)
      call p_bcast_trace_gas_type(s, pio)
      BCAST_P(s%n__num_tracer_channels )
      BCAST_P(s%n__tracer_channels     )
      BCAST_P(s%n__num_control_channels)
      BCAST_P(s%n__control_channels    )
      BCAST_P(s%n__num_flagged_channels)
      BCAST_P(s%n__flagged_channels    )
      BCAST_P(s%r__d_obs_threshold     )
      BCAST_P(s%r__d_dep_threshold     )
    end do

  end subroutine p_bcast_trace_gas

#undef  VECTOR
#undef  DERIVED
#define DERIVED type(Cloud_Detect_Type)
#define p_bcast_DERIVED p_bcast_cloud_type
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED

#undef  VECTOR
#undef  DERIVED
#define DERIVED type(Aerosol_Detect_Type)
#define p_bcast_DERIVED p_bcast_aerosol_type
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED

#undef  VECTOR
#undef  DERIVED
#define DERIVED type(Land_Sensitivity_Type)
#define p_bcast_DERIVED p_bcast_land_sens_type
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED

#undef  VECTOR
#undef  DERIVED
#define DERIVED type(Trace_Gas_Detect_Type)
#define p_bcast_DERIVED p_bcast_trace_gas_type
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_bcast_DERIVED

end module mo_cads31
