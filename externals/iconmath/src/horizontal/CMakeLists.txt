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
  iconmath-horizontal
  mo_lib_divrot.F90
  mo_lib_laplace.F90
  mo_lib_gradients.F90)

add_library(${PROJECT_NAME}::horizontal ALIAS iconmath-horizontal)

set(Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(
  iconmath-horizontal
  PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}"
  EXPORT_NAME ${PROJECT_NAME}::horizontal)

if(IM_ENABLE_LOOP_EXCHANGE)
  target_compile_definitions(iconmath-horizontal PRIVATE __LOOP_EXCHANGE)
endif()

if(IM_ENABLE_MIXED_PRECISION)
  target_compile_definitions(iconmath-horizontal PRIVATE __MIXED_PRECISION)
endif()

if(IM_ENABLE_DIM_SWAP)
  target_compile_definitions(iconmath-horizontal PRIVATE __SWAPDIM)
endif()

if(IM_ENABLE_OPENACC)
  # If _OPENACC is defined, assume that the required compiler flags are already
  # provided, e.g. in CMAKE_Fortran_FLAGS:
  if(NOT HAS_OPENACC_MACRO)
    target_compile_options(iconmath-horizontal PRIVATE ${OpenACC_Fortran_OPTIONS})
    # This make sures that unit tests (FortUTF) compiles without the need of
    # passing OpenACC compile option.
    target_link_libraries(iconmath-horizontal PRIVATE OpenACC::OpenACC_Fortran)
  endif()
endif()

target_include_directories(
  iconmath-horizontal
  PUBLIC
    # Path to the Fortran modules:
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${Fortran_MODULE_DIRECTORY}>>
    $<INSTALL_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>>
)

target_link_libraries(iconmath-horizontal PRIVATE fortran-support::fortran-support)
target_link_libraries(iconmath-horizontal PRIVATE iconmath-support)
target_link_libraries(iconmath-horizontal PRIVATE iconmath-interpolation)

install(TARGETS iconmath-horizontal EXPORT "${PROJECT_NAME}-targets")

install(DIRECTORY "${Fortran_MODULE_DIRECTORY}/" TYPE INCLUDE)
