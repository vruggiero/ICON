# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

add_library(
  iconmath-interpolation
  mo_lib_interpolation_scalar.F90
  mo_lib_interpolation_vector.F90
  mo_lib_intp_rbf.F90)

add_library(${PROJECT_NAME}::interpolation ALIAS iconmath-interpolation)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(
  iconmath-interpolation
  PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
  EXPORT_NAME ${PROJECT_NAME}::interpolation)

if(IM_ENABLE_LOOP_EXCHANGE)
  target_compile_definitions(iconmath-interpolation PRIVATE __LOOP_EXCHANGE)
endif()

if(IM_ENABLE_MIXED_PRECISION)
  target_compile_definitions(iconmath-interpolation PRIVATE __MIXED_PRECISION)
endif()

if(IM_ENABLE_DIM_SWAP)
  target_compile_definitions(iconmath-interpolation PRIVATE __SWAPDIM)
endif()

if(IM_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    target_compile_options(iconmath-interpolation
                           PRIVATE ${OpenACC_Fortran_OPTIONS})
    # This make sures that unit tests (FortUTF) compiles without the need of
    # passing OpenACC compile option.
    target_link_libraries(iconmath-interpolation PRIVATE OpenACC::OpenACC_Fortran)
  endif()
endif()

target_include_directories(
  iconmath-interpolation
  PUBLIC
    # Path to the Fortran modules:
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>>
    # Path to internal include directory
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${PROJECT_SOURCE_DIR}/include>>
)

target_link_libraries(iconmath-interpolation PRIVATE fortran-support::fortran-support)
target_link_libraries(iconmath-interpolation PRIVATE iconmath-support)

install(TARGETS iconmath-interpolation EXPORT "${PROJECT_NAME}-targets")

install(DIRECTORY "${Fortran_MODULE_DIRECTORY}/" TYPE INCLUDE)
