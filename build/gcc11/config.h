


#ifndef DUNE_COMMON_CONFIG_HH
#define DUNE_COMMON_CONFIG_HH

/* Define to 1 if you have module dune-common available */
#ifndef HAVE_DUNE_COMMON
#define HAVE_DUNE_COMMON 1
#endif



/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.11-git"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 11

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* does the standard library provide experimental::is_detected ? */
#define DUNE_HAVE_CXX_EXPERIMENTAL_IS_DETECTED 1

/* does the language support lambdas in unevaluated contexts ? */
#define DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA 1

/* does the standard library provide identity ? */
#define DUNE_HAVE_CXX_STD_IDENTITY 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the Threading Building Blocks (TBB) library */
/* #undef HAVE_TBB */




/* old feature support macros which were tested until 2.10, kept around for one more release */
/* none for 2.10 */

/* Define to ENABLE_UMFPACK if the UMFPack library is available. */
/// \deprecated Use HAVE_SUITESPARSE_UMFPACK instead
#define HAVE_UMFPACK HAVE_SUITESPARSE_UMFPACK

/* Used to call lapack functions */
#define LAPACK_NEEDS_UNDERLINE

/* If enabled certain Python modules will be precompiled */
/* #undef DUNE_ENABLE_PYTHONMODULE_PRECOMPILE */






#endif // DUNE_COMMON_CONFIG_HH



#ifndef DUNE_GEOMETRY_CONFIG_HH
#define DUNE_GEOMETRY_CONFIG_HH

/* Define to 1 if you have module dune-geometry available */
#ifndef HAVE_DUNE_GEOMETRY
#define HAVE_DUNE_GEOMETRY 1
#endif




/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "2.11-git"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR 11

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION 0





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_GEOMETRY_CONFIG_HH



#ifndef DUNE_LOCALFUNCTIONS_CONFIG_HH
#define DUNE_LOCALFUNCTIONS_CONFIG_HH

/* Define to 1 if you have module dune-localfunctions available */
#ifndef HAVE_DUNE_LOCALFUNCTIONS
#define HAVE_DUNE_LOCALFUNCTIONS 1
#endif




/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "2.11-git"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR 11

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION 0





#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif


#endif // DUNE_LOCALFUNCTIONS_CONFIG_HH



#ifndef DUNE_UGGRID_CONFIG_HH
#define DUNE_UGGRID_CONFIG_HH

/* Define to 1 if you have module dune-uggrid available */
#ifndef HAVE_DUNE_UGGRID
#define HAVE_DUNE_UGGRID 1
#endif


/* Define to the version of dune-common */
#define DUNE_UGGRID_VERSION "2.11-git"

/* Define to the major version of dune-common */
#define DUNE_UGGRID_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_UGGRID_VERSION_MINOR 11

/* Define to the revision of dune-common */
#define DUNE_UGGRID_VERSION_REVISION 0

/* begin private section */

/* see parallel/ddd/dddi.h */
/* #undef DDD_MAX_PROCBITS_IN_GID */

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
/* #undef TIME_WITH_SYS_TIME */

/* Define to 1 if UGGrid should use the complete set of green refinement rules for tetrahedra */
/* #undef DUNE_UGGRID_TET_RULESET */

/* end private section */





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_UGGRID_CONFIG_HH



#ifndef DUNE_GRID_CONFIG_HH
#define DUNE_GRID_CONFIG_HH

/* Define to 1 if you have module dune-grid available */
#ifndef HAVE_DUNE_GRID
#define HAVE_DUNE_GRID 1
#endif




/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "2.11-git"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR 2

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR 11

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION 0

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
/* #undef DUNE_ALBERTA_VERSION */

/* Define to 1 if you have mkstemp function */
#define HAVE_MKSTEMP 1







#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif

#if __has_include(<dune-uggrid-config.hh>)
  #include <dune-uggrid-config.hh>
#endif


#endif // DUNE_GRID_CONFIG_HH



#ifndef DUNE_ISTL_CONFIG_HH
#define DUNE_ISTL_CONFIG_HH

/* Define to 1 if you have module dune-istl available */
#ifndef HAVE_DUNE_ISTL
#define HAVE_DUNE_ISTL 1
#endif





/* Define to the integer type that SuperLU was compiled for
   See e.g. what int_t is defined to in slu_sdefs.h */
#define SUPERLU_INT_TYPE int

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "2.11-git"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR 2

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR 11

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION 0

/* Enable/Disable the backwards compatibility of the category enum/method in dune-istl solvers, preconditioner, etc. */
#define DUNE_ISTL_SUPPORT_OLD_CATEGORY_INTERFACE 1





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_ISTL_CONFIG_HH



#ifndef DUNE_TYPETREE_CONFIG_HH
#define DUNE_TYPETREE_CONFIG_HH

/* Define to 1 if you have module dune-typetree available */
#ifndef HAVE_DUNE_TYPETREE
#define HAVE_DUNE_TYPETREE 1
#endif




/* Define to the version of dune-typetree */
#define DUNE_TYPETREE_VERSION "2.11-git"

/* Define to the major version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MAJOR 2

/* Define to the minor version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MINOR 11

/* Define to the revision of dune-typetree */
#define DUNE_TYPETREE_VERSION_REVISION 0





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_TYPETREE_CONFIG_HH



#ifndef DUNE_FUNCTIONS_CONFIG_HH
#define DUNE_FUNCTIONS_CONFIG_HH

/* Define to 1 if you have module dune-functions available */
#ifndef HAVE_DUNE_FUNCTIONS
#define HAVE_DUNE_FUNCTIONS 1
#endif





/* Define to the version of dune-functions */
#define DUNE_FUNCTIONS_VERSION "2.11-git"

/* Define to the major version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MINOR 11

/* Define to the revision of dune-functions */
#define DUNE_FUNCTIONS_VERSION_REVISION 0





#if __has_include(<dune-localfunctions-config.hh>)
  #include <dune-localfunctions-config.hh>
#endif

#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-istl-config.hh>)
  #include <dune-istl-config.hh>
#endif

#if __has_include(<dune-typetree-config.hh>)
  #include <dune-typetree-config.hh>
#endif

#if __has_include(<dune-uggrid-config.hh>)
  #include <dune-uggrid-config.hh>
#endif


#endif // DUNE_FUNCTIONS_CONFIG_HH



#ifndef DUNE_FOAMGRID_CONFIG_HH
#define DUNE_FOAMGRID_CONFIG_HH

/* Define to 1 if you have module dune-foamgrid available */
#ifndef HAVE_DUNE_FOAMGRID
#define HAVE_DUNE_FOAMGRID 1
#endif




#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif

#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif

#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif


#endif // DUNE_FOAMGRID_CONFIG_HH



#ifndef DUNE_ALUGRID_CONFIG_HH
#define DUNE_ALUGRID_CONFIG_HH

/* Define to 1 if you have module dune-alugrid available */
#ifndef HAVE_DUNE_ALUGRID
#define HAVE_DUNE_ALUGRID 1
#endif





#define DUNE_ALUGRID_VERSION "2.11-git"

/* Define to the major version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MAJOR 2

/* Define to the minor version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MINOR 11

/* Define to the revision of dune-alugrid*/
#define DUNE_ALUGRID_VERSION_REVISION 0

/* Define to build more .cc into library */
/* #undef DUNE_ALUGRID_COMPILE_BINDINGS_IN_LIB */

/* Define if we have dlmalloc */
/* #undef HAVE_DLMALLOC */

/* Define if we have zoltan */
/* #undef HAVE_ZOLTAN */

/* Define if we have ZLIB */
#define HAVE_ZLIB 1

/* Include source file for dlmalloc */
/* #undef DLMALLOC_SOURCE_INCLUDE */

/* Define if we have thread local storage */
/* #undef HAVE_PTHREAD_TLS */

/* Define if we have pthreads */
#define HAVE_PTHREAD 1

/* Define if testgrids.hh from dune-grid have been found in docs/grids/gridfactory */
/* #undef HAVE_DUNE_GRID_TESTGRIDS */

/* Define if METIS is enabled for ALUGrid */
#define ALUGRID_HAVE_METIS 0

/* Grid type magic for DGF parser */
 
/* ALUGRID_CONFORM not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ALUGRID_CUBE not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ALUGRID_SIMPLEX not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ALUGRID_CONFORM_NOCOMM not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ALUGRID_CUBE_NOCOMM not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ALUGRID_SIMPLEX_NOCOMM not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */




#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif


#endif // DUNE_ALUGRID_CONFIG_HH



#ifndef DUNE_VTK_CONFIG_HH
#define DUNE_VTK_CONFIG_HH

/* Define to 1 if you have module dune-vtk available */
#ifndef HAVE_DUNE_VTK
#define HAVE_DUNE_VTK 1
#endif





/* Define to the version of dune-vtk */
#define DUNE_VTK_VERSION "2.11-git"

/* Define to the major version of dune-vtk */
#define DUNE_VTK_VERSION_MAJOR 2

/* Define to the minor version of dune-vtk */
#define DUNE_VTK_VERSION_MINOR 11

/* Define to the revision of dune-vtk */
#define DUNE_VTK_VERSION_REVISION 0

/* Define if you have the ZLIB library.  */
#define HAVE_VTK_ZLIB ENABLE_VTK_ZLIB





#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-localfunctions-config.hh>)
  #include <dune-localfunctions-config.hh>
#endif

#if __has_include(<dune-functions-config.hh>)
  #include <dune-functions-config.hh>
#endif

#if __has_include(<dune-spgrid-config.hh>)
  #include <dune-spgrid-config.hh>
#endif

#if __has_include(<dune-polygongrid-config.hh>)
  #include <dune-polygongrid-config.hh>
#endif

#if __has_include(<dune-alugrid-config.hh>)
  #include <dune-alugrid-config.hh>
#endif

#if __has_include(<dune-uggrid-config.hh>)
  #include <dune-uggrid-config.hh>
#endif

#if __has_include(<dune-foamgrid-config.hh>)
  #include <dune-foamgrid-config.hh>
#endif


#endif // DUNE_VTK_CONFIG_HH



#ifndef DUNE_CURVEDGEOMETRY_CONFIG_HH
#define DUNE_CURVEDGEOMETRY_CONFIG_HH

/* Define to 1 if you have module dune-curvedgeometry available */
#ifndef HAVE_DUNE_CURVEDGEOMETRY
#define HAVE_DUNE_CURVEDGEOMETRY 1
#endif





/* Define to the version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION "2.11-git"

/* Define to the major version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_MINOR 11

/* Define to the revision of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_REVISION 0





#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif

#if __has_include(<dune-localfunctions-config.hh>)
  #include <dune-localfunctions-config.hh>
#endif

#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-vtk-config.hh>)
  #include <dune-vtk-config.hh>
#endif


#endif // DUNE_CURVEDGEOMETRY_CONFIG_HH



#ifndef DUNE_CURVEDGRID_CONFIG_HH
#define DUNE_CURVEDGRID_CONFIG_HH

/* Define to 1 if you have module dune-curvedgrid available */
#ifndef HAVE_DUNE_CURVEDGRID
#define HAVE_DUNE_CURVEDGRID 1
#endif





/* Define to the version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION "2.11-git"

/* Define to the major version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_MAJOR 2

/* Define to the minor version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_MINOR 11

/* Define to the revision of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_REVISION 0





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif

#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif

#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-curvedgeometry-config.hh>)
  #include <dune-curvedgeometry-config.hh>
#endif

#if __has_include(<dune-localfunctions-config.hh>)
  #include <dune-localfunctions-config.hh>
#endif

#if __has_include(<dune-functions-config.hh>)
  #include <dune-functions-config.hh>
#endif

#if __has_include(<dune-foamgrid-config.hh>)
  #include <dune-foamgrid-config.hh>
#endif

#if __has_include(<dune-alugrid-config.hh>)
  #include <dune-alugrid-config.hh>
#endif


#endif // DUNE_CURVEDGRID_CONFIG_HH



#ifndef DUNE_MESHDIST_CONFIG_HH
#define DUNE_MESHDIST_CONFIG_HH

/* Define to 1 if you have module dune-meshdist available */
#ifndef HAVE_DUNE_MESHDIST
#define HAVE_DUNE_MESHDIST 1
#endif

/* begin hausdorff-distance
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/



/* Define to the version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION "0.1"

/* Define to the major version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_MAJOR 0

/* Define to the minor version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_MINOR 1

/* Define to the revision of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_REVISION 0

/* end hausdorff-distance
   Everything below here will be overwritten
*/



#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-functions-config.hh>)
  #include <dune-functions-config.hh>
#endif

#if __has_include(<dune-foamgrid-config.hh>)
  #include <dune-foamgrid-config.hh>
#endif

#if __has_include(<dune-vtk-config.hh>)
  #include <dune-vtk-config.hh>
#endif


#endif // DUNE_MESHDIST_CONFIG_HH



#ifndef DUNE_GMSH4_CONFIG_HH
#define DUNE_GMSH4_CONFIG_HH

/* Define to 1 if you have module dune-gmsh4 available */
#ifndef HAVE_DUNE_GMSH4
#define HAVE_DUNE_GMSH4 1
#endif





/* Define to the version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION "2.11-git"

/* Define to the major version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_MAJOR 2

/* Define to the minor version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_MINOR 11

/* Define to the revision of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_REVISION 0





#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif

#if __has_include(<dune-alugrid-config.hh>)
  #include <dune-alugrid-config.hh>
#endif

#if __has_include(<dune-foamgrid-config.hh>)
  #include <dune-foamgrid-config.hh>
#endif

#if __has_include(<dune-vtk-config.hh>)
  #include <dune-vtk-config.hh>
#endif

#if __has_include(<dune-curvedgrid-config.hh>)
  #include <dune-curvedgrid-config.hh>
#endif


#endif // DUNE_GMSH4_CONFIG_HH



#ifndef HIGHER_ORDER_BGN_CONFIG_HH
#define HIGHER_ORDER_BGN_CONFIG_HH

/* Define to 1 if you have module higher-order-bgn available */
#ifndef HAVE_HIGHER_ORDER_BGN
#define HAVE_HIGHER_ORDER_BGN 1
#endif





/* Define to the version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION "0.1"

/* Define to the major version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_MAJOR 0

/* Define to the minor version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_MINOR 1

/* Define to the revision of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_REVISION 0





#if __has_include(<dune-curvedgrid-config.hh>)
  #include <dune-curvedgrid-config.hh>
#endif

#if __has_include(<dune-meshdist-config.hh>)
  #include <dune-meshdist-config.hh>
#endif

#if __has_include(<dune-gmsh4-config.hh>)
  #include <dune-gmsh4-config.hh>
#endif

#if __has_include(<dune-vtk-config.hh>)
  #include <dune-vtk-config.hh>
#endif

#if __has_include(<dune-alugrid-config.hh>)
  #include <dune-alugrid-config.hh>
#endif

#if __has_include(<dune-foamgrid-config.hh>)
  #include <dune-foamgrid-config.hh>
#endif


#endif // HIGHER_ORDER_BGN_CONFIG_HH


#ifndef HIGHER_ORDER_BGN_CONFIG_PRIVATE_HH
#define HIGHER_ORDER_BGN_CONFIG_PRIVATE_HH


/* Name of package */
#define PACKAGE "higher-order-bgn"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "gh-zhang19@mails.tsinghua.edu.cn"

/* Define to the full name of this package. */
#define PACKAGE_NAME "higher-order-bgn"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "higher-order-bgn 0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "higher-order-bgn"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.1"



#include <higher-order-bgn-config.hh>

#endif // HIGHER_ORDER_BGN_CONFIG_PRIVATE_HH



#ifndef DUNE_COMMON_CONFIG_BOTTOM_HH
#define DUNE_COMMON_CONFIG_BOTTOM_HH



#endif // DUNE_COMMON_CONFIG_BOTTOM_HH



#ifndef DUNE_GEOMETRY_CONFIG_BOTTOM_HH
#define DUNE_GEOMETRY_CONFIG_BOTTOM_HH



#endif // DUNE_GEOMETRY_CONFIG_BOTTOM_HH



#ifndef DUNE_LOCALFUNCTIONS_CONFIG_BOTTOM_HH
#define DUNE_LOCALFUNCTIONS_CONFIG_BOTTOM_HH



#endif // DUNE_LOCALFUNCTIONS_CONFIG_BOTTOM_HH



#ifndef DUNE_UGGRID_CONFIG_BOTTOM_HH
#define DUNE_UGGRID_CONFIG_BOTTOM_HH



#endif // DUNE_UGGRID_CONFIG_BOTTOM_HH



#ifndef DUNE_GRID_CONFIG_BOTTOM_HH
#define DUNE_GRID_CONFIG_BOTTOM_HH



/* Grid type magic for DGF parser */

/* ALBERTAGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ONEDGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* YASPGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */



#endif // DUNE_GRID_CONFIG_BOTTOM_HH



#ifndef DUNE_ISTL_CONFIG_BOTTOM_HH
#define DUNE_ISTL_CONFIG_BOTTOM_HH



#endif // DUNE_ISTL_CONFIG_BOTTOM_HH



#ifndef DUNE_TYPETREE_CONFIG_BOTTOM_HH
#define DUNE_TYPETREE_CONFIG_BOTTOM_HH



#endif // DUNE_TYPETREE_CONFIG_BOTTOM_HH



#ifndef DUNE_FUNCTIONS_CONFIG_BOTTOM_HH
#define DUNE_FUNCTIONS_CONFIG_BOTTOM_HH



#endif // DUNE_FUNCTIONS_CONFIG_BOTTOM_HH



#ifndef DUNE_FOAMGRID_CONFIG_BOTTOM_HH
#define DUNE_FOAMGRID_CONFIG_BOTTOM_HH



#endif // DUNE_FOAMGRID_CONFIG_BOTTOM_HH



#ifndef DUNE_ALUGRID_CONFIG_BOTTOM_HH
#define DUNE_ALUGRID_CONFIG_BOTTOM_HH



#endif // DUNE_ALUGRID_CONFIG_BOTTOM_HH



#ifndef DUNE_VTK_CONFIG_BOTTOM_HH
#define DUNE_VTK_CONFIG_BOTTOM_HH



#endif // DUNE_VTK_CONFIG_BOTTOM_HH



#ifndef DUNE_CURVEDGEOMETRY_CONFIG_BOTTOM_HH
#define DUNE_CURVEDGEOMETRY_CONFIG_BOTTOM_HH



#endif // DUNE_CURVEDGEOMETRY_CONFIG_BOTTOM_HH



#ifndef DUNE_CURVEDGRID_CONFIG_BOTTOM_HH
#define DUNE_CURVEDGRID_CONFIG_BOTTOM_HH



#endif // DUNE_CURVEDGRID_CONFIG_BOTTOM_HH



#ifndef DUNE_MESHDIST_CONFIG_BOTTOM_HH
#define DUNE_MESHDIST_CONFIG_BOTTOM_HH



#endif // DUNE_MESHDIST_CONFIG_BOTTOM_HH



#ifndef DUNE_GMSH4_CONFIG_BOTTOM_HH
#define DUNE_GMSH4_CONFIG_BOTTOM_HH



#endif // DUNE_GMSH4_CONFIG_BOTTOM_HH



#ifndef HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH
#define HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH



#endif // HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH

