!
!+ Interface to the RADAR observation operator
!
module mo_radar_obs
!
! Description:
!   Interface to radar observation operator.
!   It organizes the reading and processing of
!   RADAR observation reports. It calls the 
!   specific routines of the RADAR observation 
!   operator which are compiled conditionally.
!
! Method:
!    This module contains the following module procedures:
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2008.
! Software Standards:
!
!------------------------------------------------------------------------------
  use mo_kind,         only: wp                         ! working precision kind parameter
  use mo_namelist,     only: nnml,                     &! Fortran unit used for namelist file
                             position_nml,             &! subroutine: position namelist file
                             POSITIONED, MISSING        ! return values from position_nml

  use mo_mpi_dace,     only: dace,                     &! MPI group info
                             p_bcast                    ! generic broadcast routine

  use mo_exception,    only: finish                     ! abort
  use mo_atm_grid,     only: t_grid                     ! grid derived type
  use mo_atm_state,    only: t_atm                      ! atmospheric state data-type

#ifdef USE_RADAR_OBS
  use radar_interface, only: iunit_nml,                &! namelist file Fortran unit
                             p_grid,                   &! pointer to t_grid
                             p_atm,                    &! pointer to t_atm
                             get_model_config_for_radar,&! set communication parameters
                             get_grid_variables,       &
                             get_atm_variables
  use data_radar_mie,  only: itype_gscp                 ! type of grid-scale precipitation physics
  use src_radar,       only: organize_radar,           &! Interface to radar operator
                             read_obs_rad               ! read the radar data in DWD NetCDF-Format
  use data_radar,      only: ireals,                   &! real (dp) used in EMVORADO
                             nradsta,                  &! number of radar stations
                             num_compute,              &! number of compute PEs
                             radar_data_type,          &! array of radar data (local PE) 
                             rs_meta,                  &! radar meta data for each station
                             rs_data                    ! radar data
#endif
  implicit none
!=============================================================================

  !----------------
  ! Public entities
  !----------------
  private
  public :: init_radar_obs     ! initialize RADAR observation operator modules
  public :: read_radar_obs     ! organize reading and distribution of reports
  public :: radar2dace         ! fill radar data into DACE observation type
!=============================================================================
contains
!=============================================================================

  subroutine init_radar_obs (grid)
    !type(t_grid) ,intent(in) :: grid
    type(t_grid) ,target ,intent(in) :: grid
    integer                          :: ierr
    
    integer :: nnew = 3          ! Aendern!!! Nur fuer technischen Test! [EB]
                                 ! in cosmo_refl_offline: nnew = 1         ! nnew has to be 1 ???

    !----------------------------------------------
    ! initialisation for the radar forward operator
    !----------------------------------------------
#ifndef USE_RADAR_OBS
    call finish('init_radar_obs','radar operator sources not linked')
#else
    !--------------------------------------
    ! provide unit number for namelist file
    !--------------------------------------
    iunit_nml =  nnml
    p_grid    => grid
    !---------------------------------------------------------------------
    ! Initialize the environment for the radar forward operator EMRADSCOPE
    !---------------------------------------------------------------------
    CALL get_model_config_for_radar (grid)
    !---------------------------------------------------------------------
    ! Read EMVORADO offline specific namelist REFL_OFFLINE
    !---------------------------------------------------------------------
    call read_nml_refl_offline
    !---------------------------------------------------------------------
    ! Connect grid variables with associated pointers in forward operator
    !---------------------------------------------------------------------
    call get_grid_variables (p_grid)
    !---------------------------------------------------------------------
    ! Initialization of radar forward operator, reading the namelist
    !---------------------------------------------------------------------
    if (dace% lpio) call position_nml ('RADARSIM_PARAMS', status=ierr)
    call organize_radar ('init', nnew)
#endif
  end subroutine init_radar_obs

!=============================================================================

#ifdef USE_RADAR_OBS
  subroutine read_nml_refl_offline ()
    
    integer       :: itype_gscp_d  ! type of grid-scale precipitation physics, default 3
    !integer       :: itype_gscp    ! type of grid-scale precipitation physics
    integer       :: ierr

    namelist /REFL_OFFLINE/ itype_gscp

    !-------------------------------------------------------------------------------
    !- Section 1: Initialize the default variables
    !-------------------------------------------------------------------------------
    itype_gscp_d = 3

    !-------------------------------------------------------------------------------
    !- Section 2: Initialize variables with defaults
    !-------------------------------------------------------------------------------
    itype_gscp = itype_gscp_d

    !-------------------------------------------------------------------------------
    !- Section 3: Read namelist variables and replace default variables if necessary
    !-------------------------------------------------------------------------------
    !----------------------------
    ! read on processor p_io only
    !----------------------------
    if (dace% lpio) then
      !----------------------------------
      ! search for namelist group in file
      !----------------------------------
      call position_nml ('REFL_OFFLINE', status=ierr)
      select case (ierr)
      case (POSITIONED)
        read (nnml ,nml=REFL_OFFLINE)
      case (MISSING)
      case default
        call finish('read_nml_refl_offline','ERROR reading namelist group /REFL_OFFLINE/')
      end select
    end if
    !---------------------------
    ! broadcast namelist entries
    !---------------------------
    call p_bcast (itype_gscp, dace% pio)
    !--------------------------
    ! printout namelist entries
    !--------------------------
    if (dace% lpio) then
      print *, 'namelist /REFL_OFFLINE/: '
      print *, ' '
      print *, '  itype_gscp = ', itype_gscp
    end if   

  end subroutine read_nml_refl_offline
#endif

!=============================================================================

  subroutine read_radar_obs ()
    integer :: nnew = 3          ! Aendern!!! Nur fuer technischen Test! [EB]
                                 ! in cosmo_refl_offline: nnew = 1         ! nnew has to be 1 ???
    integer                          :: ista, ista_par    ! loop variables over radar stations
    integer                          :: p_pe_rad          ! PE that reads the radar observations
!    type(radar_data_type)            :: rs_obs(nradsta)   ! type to store radar observations
#ifdef USE_RADAR_OBS
    real (kind=ireals)               :: time_mod          ! seconds since model start
#endif

#ifndef USE_RADAR_OBS
    call finish('read_radar_obs','radar operator sources not linked')
#else   

    if (dace% lpio)    print *, "test: read_radar_obs subroutine"
    if (dace% lpio)    print *, "nradsta = ", nradsta

    !---------------------------------------------------------------
    ! Read radar data from station number "ista" at time "time_mod":
    ! Data will be stored in rs_data of type radar_data_type
    !---------------------------------------------------------------

    time_mod = 0
    ! loop over num_compute-sized chunks of radar stations:
    ! (see chunkloop and istaloop in src_radar.f90)
    do ista_par = 1, nradsta, num_compute   
      do ista = ista_par, min(ista_par+num_compute-1, nradsta)
        ! Define which pe should read the radar data 
        p_pe_rad = mod(ista-1,num_compute)

        ! Read radar data for current time from cdfin file
        if (dace% pe == p_pe_rad) call read_obs_rad(ista, time_mod, 'read')
        ! Apply superobbing
        
      enddo
    end do

#endif

  end subroutine read_radar_obs

  subroutine compute_radar_mod (atm)
    type(t_atm)  ,target ,intent(in) :: atm
!    type(t_grid) ,target ,intent(in) :: grid   ! ueberhaupt noetig??? [EB]

#ifndef USE_RADAR_OBS
    call finish('compute_radar_mod','radar operator sources not linked')
#else 
    p_atm => atm

    call get_atm_variables (p_atm)
!     call organize_radar ('compute', nnew)
#endif

  end subroutine compute_radar_mod

!=============================================================================

  subroutine radar2dace
  !-------------------------------------------------
  ! fill radar data into DACE observation type t_obs
  !  +++ work in progress +++
  !-------------------------------------------------

    integer :: ista, ista_par    ! loop variables over radar stations
    integer :: p_pe_rad          ! processor for this radar station
    integer :: it, iel, iaz, ira ! time,elevation,azimuth,range index
    integer :: nt, nel, naz, nra ! number of: times,elevations,azimuths,ranges
    integer :: ind               ! index in adwind_obs array
    integer :: n_obs             ! observation counter

#ifndef USE_RADAR_OBS
    call finish('radar2dace','radar operator sources not linked')
#else
    !--------------------------------------------------
    ! loop over radar stations, select processor to use
    !--------------------------------------------------
    do ista_par = 1, nradsta, dace% npe
      do ista = ista_par, min(ista_par+dace% npe-1, nradsta)
        p_pe_rad = mod(ista-1, dace% npe)
        if (dace% pe == p_pe_rad) then
        !--------------------------------------
        ! fetch some meta data for this station
        !--------------------------------------
        nra = rs_meta (ista)% nra  ! number of ranges
        naz = rs_meta (ista)% naz  ! number of azimuths
        nel = rs_meta (ista)% nel  ! number of elevations
        !-----------------------------------
        ! loop over time slots
        ! +++ currently only time slot 1 +++
        !-----------------------------------
        do it = 1, rs_meta (ista)% nobs_times
!print *,dace% pe,'### ista it obs_times',ista, it, rs_meta (ista)% obs_times(it)
          if (rs_meta (ista)% obs_times(it) == 0) then
            !---------------------
            ! loop over elevations
            !---------------------
            do iel = 1, nel
!print *,dace% pe,'### ista it obs_times iel el',ista, it, rs_meta (ista)% obs_times(it), &
!                                                iel, rs_meta (ista)% el_arr(iel)
              n_obs = 0
              !------------------------------
              ! loop over azimuths and ranges
              !------------------------------
              do iaz = 1, naz
                do ira = 1, nra
                  ind = iaz + (ira-1) * naz + (iel-1) * (naz * nra)
if (iel == 1 .and. iaz == 1) then
if (rs_data(ista) % radwind_obs (ind) < 1.e35_wp) &
print *,dace% pe,'### data', rs_data(ista) % radwind_obs (ind)
endif
                end do
              end do
            end do
          endif
        end do
        endif
      end do
    end do
#endif
  end subroutine radar2dace

!==============================================================================
end module mo_radar_obs
