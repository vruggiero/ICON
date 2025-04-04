#!/bin/ksh
#
# gen_coast.com -  Script that generates coastal area from a given land sea mask
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
set -ex
#
# *********** Data for Program *************************
# *** This needs to be EDITED  *************************
DATA=/mnt/lustre01/work/gg0302/g260122/HD/forcing  # Data directory
DNAM=land_cd3.nc                                  # Land sea mask file
dscr=/mnt/lustre01/pf/g/g260122/hdmodel/util       # Script/Program directory
DMAX=50               # Size of coastal zone in [km]
# ******************************************************
#
#  COMPILER settings

LPATH="-L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib"
LIBS="-lnetcdf"

NETCDF_LIB="`/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --flibs`"
NC_INCLUDE="-I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include"
#
# *** ifort klappt nicht, ebenso mpiifort
F90="mpif90 -c -fpp ${NC_INCLUDE}"

FFLAGS="-Os -no-vec ${NC_INCLUDE}"
LDPAR="mpif90"
LDFLG="-lz -shared-intel ${LPATH} ${LIBS} ${NETCDF_LIB} ${NC_INCLUDE}"
#
# *** set up modules for Bull compiler
source /sw/rhel6-x64/etc/profile.mistral
  [[ -n `whence ifort` ]] && module unload intel
   module load intel/17.0.0
  [[ -n `whence mpif90` ]] && module unload bullxmpi_mlx
   module load mxm/3.3.3002
   module load fca/2.5.2393
   module load bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
#
# **************** Compile **************************
cd $TMPDIR
${F90} ${FFLAGS} $dscr/mo_interpol.f90 $dscr/generate_coast.f90
${LDPAR} ${LDFLG} -o gen_coast.exe mo_interpol.o generate_coast.o
#
cdo setvar,mask ${DATA}/$DNAM lsm.nc
./gen_coast.exe $DMAX

#

