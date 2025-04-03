! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Main program for the ICON atmospheric model

MODULE mo_icon_testbed

  USE mo_exception,           ONLY: finish
  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_icon_testbed_config, ONLY: testbed_model, null_model, test_coupler_model, &
    & test_jitter_model, test_halo_communication, test_netcdf_read_model,          &
    & test_gather_communication, test_exchange_communication,                      &
    & test_bench_exchange_data_mult
  USE mo_icon_testbed_nml,    ONLY: read_icon_testbed_namelist

#ifndef __NO_ICON_ATMO__
  USE mo_test_communication,  ONLY: test_communication
  USE mo_test_jitter,         ONLY: test_jitter
  USE mo_test_netcdf_read,    ONLY: test_netcdf_read
#endif

#ifndef __NO_ICON_COMIN__
  USE mo_mpi,               ONLY: p_comm_comin
  USE comin_host_interface, ONLY: mpi_handshake_dummy
#endif
!-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: icon_testbed

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE icon_testbed(testbed_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: testbed_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_icon_testbed:icon_testbed"

    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    CALL read_icon_testbed_namelist(testbed_namelist_filename)

#ifndef __NO_ICON_COMIN__
    ! we dont participate at comin (yet) but we need to be friendly and shake hands
    CALL mpi_handshake_dummy(p_comm_comin)
#endif

    SELECT CASE(testbed_model)
    
    CASE(null_model)
      ! do nothing
      RETURN


#ifndef __NO_ICON_ATMO__

    CASE(test_halo_communication, test_gather_communication, &
         test_exchange_communication, test_bench_exchange_data_mult)
      CALL test_communication(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_jitter_model)
      CALL test_jitter(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_netcdf_read_model)
      CALL test_netcdf_read(testbed_namelist_filename,shr_namelist_filename)
#endif

    CASE default
      CALL finish(method_name, "Unrecognized testbed_model")

    END SELECT    
   

  END SUBROUTINE icon_testbed
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed

