EMVORADO sources to serve as a submodule in ICON
================================================


Maintainer:      Ulrich Blahak, DWD
Creation date:   8.4.2020
Last change:     8.4.2020
Contact:         ulrich.blahak@dwd.de

Administration: https://gitlab.dkrz.de/dwd-sw/emvorado-for-icon/

Short intro:
============

This submodule contains the necessary sources for the
efficient modular volume scan radar forward operator (EMVORADO) to
be coupled to ICON.

To configure ICON to include EMVORADO, simply add the configuration flag
  --enable-emvorado
to your config wrapper call, i.e.,
  $> ./config/dwd/xce.cray --enable-emvorado

EMVORADO is also coupled to other NWP-models (COSMO, in future maybe WRF) 
and shares most of its core modules among these models. Many parts of
its code do not depend on the specific model, but on the specifics
of radars and their volume scans. Moreover, all these NWP-models are MPI-
parallel using similar concepts of domain decomposition, and their
grids and model variables are in principle 3D.

To adress the specific differences of the grids and data structures,
at the same time preserving computational accuracy and minimizing
runtime, there are model specific interfaces, which work in both directions:
The model can call procedures provided by EMVORADO for initialization and
radar forward simulation, while EMVORADO can directly use the model
internal data structures (model fields, grid information types) for
interpolation from model grids to radar bins (lat/lon/height).

For example, the model state variables are not copied to any transfer
variables or interpolated to some sort of intermediate grid, but
just directly linked to EMVORADO-internal 3D pointers, which is very efficient.

However, this requires that the EMVORADO interface in the model code
(ICON: src/data_assimilation/interfaces/radar_interface.f90)
provides model-specific and efficient interpolation code from the model grid positions
to the radar bins. This code fills the bodies of interpolation procedures which have
a model-independent interface and can be generically called by EMVORADO routines.
In ICON, this needs access to the ICON grid data structures in the radar_interface.f90.

In summary, there are model-independend EMVORADO sources and model specific
interfaces. All code which explicitly references model code in any fashion
is contained in the specific interfaces.

In case of ICON, the model specific interface modules reside in

${ICON_ROOT_DIR}/src/data_assimilation/interfaces/ ...
       ... radar_interface.f90
       ... mo_emvorado_config.f90
       ... mo_emvorado_init.f90
       ... mo_emvorado_interface.f90

The model-independent EMVORADO sources, which will be compiled into ICON
by its build process, are in the subdirectory

./src_emvorado/

The model-dependend sources (currently only 1 file, radar_mpi_init_icon.f90)
are in subdirectory

./src_iface_icon/

For the generation of netcdf-feedback files for data assimilation, EMVORADO
relies on a few type declarations from modules of
DWD's DACE (data assimilation coding environment) code:

o If ICON is configured with the dace_icon submodule (--enable-dace), the
  corresponding modules are present and used from there.
o IF ICON is not configured that way (--enable-dace not set, or --disable-dace set),
  we use these modules from the dace_icon submodule anyways, triggered by an extra entry in
  ${ICON_ROOT_DIR}/icon.mk.in.

The directory structure of the present submodule is a subset of the original
EMVORADO repository (see below "History") and is kept that way to make bookkeeping and backmerging
to this repository easier.


Documentation:
==============

./DOC/TEX/emvorado-userguide.pdf

and the references therein. This user guide is also available online:

www.cosmo-model.org --> Model Documentation --> EMVORADO


History:
========

The present git repository

  git@gitlab.dkrz.de:dwd-sw/emvorado-for-icon.git

has been created by pulling the
former branch "submodule-code-for-icon-nwp-dev" from the
full EMVORADO repository

  git@gitlab.dkrz.de:dace_projects/emvorado-package.git

to here and merging it into master. Many for ICON unnecessary things
from the full repo have been removed, while preserving the original
directory structure. In this way we preserve the
history of the full EMVORADO repository and enable EMVORADO developers
to exchange developments between the repositories.
At the same time, we limit the amount of data in this submodule to a minimum.


To bring developments back from ICON side (for developers):
===========================================================

In ICON:
--------

cd externals/emvorado
git checkout -b dev/from-icon
git add ...
git commit
git push -u origin dev/from-icon

In emvorado-for-icon.git: merge commits from dev/from-icon to master via roundtrip cosmo-emvorado-packate.git:
--------------------------------------------------------------------------------------------------------------

git checkout -b dev/from-icon origin/dev/from-icon
git pull
git remote add cosmo-emvorado git@gitlab.dkrz.de:dace_projects/emvorado-package.git
git fetch cosmo-emvorado
git checkout -b cosmo-master cosmo-emvorado/master
git cherry-pic <commit1..commit2>  # the commits from dev/from-icon to transfer to cosmo-master
git push cosmo-emvorado HEAD:master

git checkout master
git merge cosmo-master

git remote rm cosmo-emvorado
git branch -D cosmo-master

git branch -D dev/from-icon
git push --delete origin dev/from-icon

