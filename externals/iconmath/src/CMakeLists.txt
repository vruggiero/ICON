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

add_subdirectory(support)

if(BUILD_ICONMATH_INTERPOLATION)
  add_subdirectory(interpolation)
endif()

if(BUILD_ICONMATH_HORIZONTAL)
  add_subdirectory(horizontal)
endif()
