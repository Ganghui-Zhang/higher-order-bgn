General requirements:
======================================
- C++ compiler supporting the C++17 language standard, e.g. g++ >= 7 or clang >= 5
- CMake >= 3.13
- SuiteSparse / UMFPack

Getting all Dune dependencies:
======================================
The code depends on several dune modules, publicly available
- [dune-common](https://gitlab.dune-project.org/core/dune-common.git)
- [dune-geometry](https://gitlab.dune-project.org/core/dune-geometry.git)
- [dune-grid](https://gitlab.dune-project.org/core/dune-grid.git)
- [dune-istl](https://gitlab.dune-project.org/core/dune-istl.git)
- [dune-localfunctions](https://gitlab.dune-project.org/core/dune-localfunctions.git)
- [dune-functions](https://gitlab.dune-project.org/staging/dune-functions.git)
- [dune-typetree](https://gitlab.dune-project.org/staging/dune-typetree.git)
- [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid.git)
- [dune-alugrid](https://gitlab.dune-project.org/extensions/dune-alugrid.git)
- [dune-foamgrid](https://gitlab.dune-project.org/extensions/dune-foamgrid.git)
- [dune-curvedgeometry](https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgeometry.git)
- [dune-curvedgrid](https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgrid.git)
- [dune-gmsh4](https://gitlab.mn.tu-dresden.de/iwr/dune-gmsh4.git)
- [dune-vtk](https://gitlab.mn.tu-dresden.de/iwr/dune-vtk.git)

And on the *new* dune modules introduced in the manuscript:
- [dune-meshdist](https://gitlab.dune-project.org/extentsions/dune-meshdist)

Download all modules in a common root directory, called `DUNE_TEST`, change to that directory
and run the command `dunecontrol` to configure and build all dependencies:

./dune-common/bin/dunecontrol all

Compiling the numerical tests using the gcc11.opts file
======================================
./dune-common/bin/dunecontrol --opts=gcc11.opts --only=higher-order-bgn all

Running the numerical tests
======================================
When build with cmake, all numerical tests executables are placed in the directory `build/gcc11/src`.
Thus, in order to run a numerical tests `<example>`, simply call

./build/gcc11/src/example

Explicit numerical tests (Need to adjust the local assembler matrix)
======================================
- The experiment order of convergence for mean curvature flow:
  Complie eoc_curve_test_ref.cc and eoc_test.hh 
- The experiment order of convergence for surface diffusion:
  Complie eoc_curve_test.cc and eoc_test.hh
- The evolution test for mean curvature flow and surface diffusion of curve:
  Complie higher-order-bgn_curve_cc and Curve.hh
- The evolution test for mean curvature flow and surface diffusion of surface:
  Complie higher-order-bgn_surface_cc and Surface.hh
- The structure-preserving evolution test for surface diffusion of curve:
  Complie higher-order-bgn_curve_cc and Curve_SP.hh
- The structure-preserving evolution test for surface diffusion of surface:
  Complie higher-order-bgn_surface_cc and Surface_SP.hh
