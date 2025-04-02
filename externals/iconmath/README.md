<!--
ICON

---------------------------------------------------------------
Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: CC-BY-4.0
---------------------------------------------------------------
-->

# ICONMATH library
This repository is an external library of ICON collecting modules of ICON supporting and enabling math operations. The library is divided into three components:
- `libiconmath-support`: contains general modules that are used in ICON to support math operations but are independent of the data types in ICON.
- `libiconmath-horizontal`: contains modules which compute several mathematical operators such as divergence, gradient, Laplacian etc.
- `libiconmath-interpolation`: contains modules that interpolates and reconstructs scalar and vector fields used in ICON.

[![pipeline status](https://gitlab.dkrz.de/icon-libraries/libiconmath/badges/master/pipeline.svg)](https://gitlab.dkrz.de/icon-libraries/libiconmath/pipelines/latest?ref=master)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://gitlab.dkrz.de/icon-libraries/libiconmath/-/blob/master/LICENSES/BSD-3-Clause.txt)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libiconmath/-/blob/master/LICENSES/CC-BY-4.0.txt)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libiconmath/-/blob/master/LICENSES/CC0-1.0.txt)

## Requirements
The following packages/libraries are required for `libiconmath` libraries
- Fortran compiler
- CMake 3.18+
- [libfortran-support](https://gitlab.dkrz.de/icon-libraries/libfortran-support) to access some utility functions for Fortran

The following requirements are optional

- `fprettify` for Fortran code formatting
- `OpenACC` to run the code on GPU nodes
- [REUSE](https://reuse.software) v3.0.0+ for licensing

---
# For Users

## How to install `libiconmath`?

1. Clone the repository
```bash
git clone https://gitlab.dkrz.de/icon-libraries/libiconmath.git
cd libiconmath
```

2. Create a build directory
```bash
mkdir build
cd build
```

3. Configure the library
```bash
cmake ..
```
For GPU support with OpenACC, use the following command
```bash
cmake -DIM_ENABLE_OPENACC=ON ..
```

4. Compile the library
```bash
make
```

5. Install the library
```bash
make install
```
`make install` needs to be run with sudo if the installation directory is not writable by the user. To install the library in a specific directory, use the following command
```bash
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX ..
make install
```
where `$PREFIX` is the installation directory.

6. Check the installation
```bash
ls $PREFIX/lib
```
where `$PREFIX` is the installation directory. The library files should be in the `lib` directory.

7. Check the module files
```bash
ls $PREFIX/include
```
where `$PREFIX` is the installation directory. The module files should be in the `mod` directory.

## How to use `libiconmath` in your project?

### What modules are in the `libiconmath-support` library?
The `libiconmath-support` library includes some general modules that are used in ICON to support math operations but are independent of the data types in ICON. Here is a list of the supported modules.
- `mo_gridman_constants`: constants relevant for implementational issues
- `mo_lib_grid_geometry_info`: grid-related basic structures and geometry parameters
- `mo_lib_loopindices`: determines start and end indices of do loops
- `mo_math_constants`: main mathematical constants
- `mo_math_types`: basis math types
- `mo_math_utilities`: various mathematical algorithms

### What modules are in the `libiconmath-interpolation` library?
The `libiconmath-interpolation` library includes modules that interpolates and reconstructs scalar and vector fields used in ICON. Here is a list of the supported modules.
- `mo_lib_interpolation_scalar`: contains the implementation of interpolation and reconstruction for scalar fields
- `mo_lib_interpolation_vector`: contains the implementation of interpolation and reconstruction for vector fields
- `mo_lib_intp_rbf`: perform vector RBF reconstruction

### What modules are in the `libiconmath-horizontal` library?
The `libiconmath-horizontal` library includes modules which compute several mathematical operators such as divergence, gradient, Laplacian etc.
- `mo_lib_divrot`: contains the implementation of divergence, rotation and linear reconstruction
- `mo_lib_gradient`: contains the implementation of mathematical gradient operators
- `mo_lib_laplace`: contains the implementation of computing Laplacian of scalar and vector fields

### How to link modules from math libraries (same for `libiconmath-support`, `libiconmath-interpolation`, and `libiconmath-horizontal`)
- The individual libraries will be built into `libiconmath/build/src/support/<library_name>.dylib` (or `.so` depending on the platform)
- The module files (Fortran) for a `<component>` of the library (i.e. `support`, `interpolation` or `horizontal`) will be built under `libiconmath/build/src/<component>/mod/`

### Include the library in your CMakeLists.txt
```cmake
find_package(iconmath REQUIRED)
target_link_libraries(your_target PRIVATE iconmath::support iconmath::interpolation iconmath::horizontal)
```

### Include the module files in your Fortran code, as an example,
```fortran
use mo_math_types, only: t_tangent_vectors
use mo_lib_interpolation_scalar, only: verts2edges_scalar_lib
use mo_lib_divrot, only: div_lib
```
---
# For developers

## Some notes for developers
- The `libiconmath` library is only configured by CMake.
  - Tips and standards on CMake [ICON developer wiki/CMake recommendations and requirements](https://gitlab.dkrz.de/icon/wiki/-/wikis/CMake-recommendations-and-requirements)
- The `libiconmath` library uses `fprettify` for formatting Fortran codes. Run `make format` before you commit.
- The `libiconmath` library is unit tested. (work in progress) All merge request changes are preferable to have a unit test. See [icon-c/Wiki/Testing and building of ICON C/Unit test frameworks](https://gitlab.dkrz.de/icon/icon-c/-/wikis/ICON-C-Phase-0/Testing-and-building-of-ICON-C#unit-test-frameworks) for more information on unit testing.
  - __Fortran__ unit tests are written in [FortUTF](https://github.com/artemis-beta/FortUTF). [Assertions Documentation](https://github.com/artemis-beta/FortUTF/blob/main/docs/assertions.md)
- Fortran preprocessing is automatically applied for files with `.F90` extensions. See [\#4](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/issues/4) for more details.
- Make sure your code is tested.
- Check license by `reuse lint`. Check [requirements](#requirements).

## How to add modules in `libiconmath-support`, `libiconmath-interpolation` or `libiconmath-horizontal`?
- For each `<component>`, the steps are the same. The only difference is the directory where the module file should be placed.
1. Put your module file under `src/support/`, `src/interpolation/` or `src/horizontal/` depending on the component.
2. Add your file to the library configuration list at `src/<component>/CMakeLists.txt`.
```
add_library(
  iconmath-<component>
  ...
# Add your module to the list. The list is in alphabetical order.
  ...
)
```
## Check linting and formatting
### Check the code format
```bash
make format
```
### Check the license
- Use the following command to check whether the license header is added correctly
```bash
reuse lint
```
- More details can be found in the [reuse tool](https://github.com/fsfe/reuse-tool) GitHub page.
- The following license should be used for the library
  - Code snippets should have license BSD-3-Clause
  - Documentations should have license CC-BY-4.0
  - Files unrelated to the library itself should have license CC0-1.0
### Check the CMake style
- Use the following command to check the CMake style
```bash
cmake-lint cmake/*.cmake **/CMakeLists.txt
```
- Use the following command to fix the CMake style issues
```bash
cmake-format -i cmake/*.cmake **/CMakeLists.txt
```
- More details can be found in the [cmake-format](https://github.com/cheshirekow/cmake_format) GitHub page.
### Check typos
- Use the following command to check typos in the code
```bash
typos
```
- More details can be found in the [typos](https://github.com/crate-ci/typos) GitHub page.
### Check OpenACC style
- Clone the repository [ICON OpenAcc Beautifier](https://gitlab.dkrz.de/dwd-sw/icon-openacc-beautifier)
```bash
git clone https://gitlab.dkrz.de/dwd-sw/icon-openacc-beautifier.git
```
- Use the following command to check the OpenACC style
```bash 
python3 icon-openacc-beautifier/main.py src/
```

## How to contribute
Please open a merge request and select one of our templates: __[feature/bugfix]__. Detailed instructions on how to proceed are provided there.

## Contact
This repository is mainly maintained by the following maintainers:
- __Pradipta Samanta__ (samanta@dkrz.de)
- __Yen-Chen Chen__ (yen-chen.chen@kit.edu)
- __Sergey Kosukhin__ (sergey.kosukhin@mpimet.mpg.de)

This repository is owned by the `icon-libraries` group, contacts about general ICON library questions:
- __Will Sawyer__ (william.sawyer@cscs.ch)
- __Florian Prill__ (florian.prill@dwd.de)
- __Luis Kornblueh__ (luis.kornblueh@mpimet.mpg.de)
