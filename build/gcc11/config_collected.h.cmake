


#ifndef DUNE_COMMON_CONFIG_HH
#define DUNE_COMMON_CONFIG_HH

/* Define to 1 if you have module dune-common available */
#ifndef HAVE_DUNE_COMMON
#cmakedefine01 HAVE_DUNE_COMMON
#endif



/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "${DUNE_COMMON_VERSION}"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR ${DUNE_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR ${DUNE_COMMON_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION ${DUNE_COMMON_VERSION_REVISION}

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL ${DUNE_MINIMAL_DEBUG_LEVEL}

/* does the standard library provide experimental::is_detected ? */
#cmakedefine DUNE_HAVE_CXX_EXPERIMENTAL_IS_DETECTED 1

/* does the language support lambdas in unevaluated contexts ? */
#cmakedefine DUNE_HAVE_CXX_UNEVALUATED_CONTEXT_LAMBDA 1

/* does the standard library provide identity ? */
#cmakedefine DUNE_HAVE_CXX_STD_IDENTITY 1

/* Define if you have a BLAS library. */
#cmakedefine HAVE_BLAS 1

/* Define if you have LAPACK library. */
#cmakedefine HAVE_LAPACK 1

/* Define to 1 if you have the Threading Building Blocks (TBB) library */
#cmakedefine HAVE_TBB 1




/* old feature support macros which were tested until 2.10, kept around for one more release */
/* none for 2.10 */

/* Define to ENABLE_UMFPACK if the UMFPack library is available. */
/// \deprecated Use HAVE_SUITESPARSE_UMFPACK instead
#define HAVE_UMFPACK HAVE_SUITESPARSE_UMFPACK

/* Used to call lapack functions */
#cmakedefine LAPACK_NEEDS_UNDERLINE

/* If enabled certain Python modules will be precompiled */
#cmakedefine DUNE_ENABLE_PYTHONMODULE_PRECOMPILE






#endif // DUNE_COMMON_CONFIG_HH



#ifndef DUNE_GEOMETRY_CONFIG_HH
#define DUNE_GEOMETRY_CONFIG_HH

/* Define to 1 if you have module dune-geometry available */
#ifndef HAVE_DUNE_GEOMETRY
#cmakedefine01 HAVE_DUNE_GEOMETRY
#endif




/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "${DUNE_GEOMETRY_VERSION}"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR ${DUNE_GEOMETRY_VERSION_MAJOR}

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR ${DUNE_GEOMETRY_VERSION_MINOR}

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION ${DUNE_GEOMETRY_VERSION_REVISION}





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_GEOMETRY_CONFIG_HH



#ifndef DUNE_LOCALFUNCTIONS_CONFIG_HH
#define DUNE_LOCALFUNCTIONS_CONFIG_HH

/* Define to 1 if you have module dune-localfunctions available */
#ifndef HAVE_DUNE_LOCALFUNCTIONS
#cmakedefine01 HAVE_DUNE_LOCALFUNCTIONS
#endif




/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "${DUNE_LOCALFUNCTIONS_VERSION}"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR ${DUNE_LOCALFUNCTIONS_VERSION_MAJOR}

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR ${DUNE_LOCALFUNCTIONS_VERSION_MINOR}

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION ${DUNE_LOCALFUNCTIONS_VERSION_REVISION}





#if __has_include(<dune-geometry-config.hh>)
  #include <dune-geometry-config.hh>
#endif


#endif // DUNE_LOCALFUNCTIONS_CONFIG_HH



#ifndef DUNE_UGGRID_CONFIG_HH
#define DUNE_UGGRID_CONFIG_HH

/* Define to 1 if you have module dune-uggrid available */
#ifndef HAVE_DUNE_UGGRID
#cmakedefine01 HAVE_DUNE_UGGRID
#endif


/* Define to the version of dune-common */
#define DUNE_UGGRID_VERSION "${DUNE_UGGRID_VERSION}"

/* Define to the major version of dune-common */
#define DUNE_UGGRID_VERSION_MAJOR ${DUNE_UGGRID_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_UGGRID_VERSION_MINOR ${DUNE_UGGRID_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_UGGRID_VERSION_REVISION ${DUNE_UGGRID_VERSION_REVISION}

/* begin private section */

/* see parallel/ddd/dddi.h */
#cmakedefine DDD_MAX_PROCBITS_IN_GID ${UG_DDD_MAX_MACROBITS}

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#cmakedefine TIME_WITH_SYS_TIME 1

/* Define to 1 if UGGrid should use the complete set of green refinement rules for tetrahedra */
#cmakedefine DUNE_UGGRID_TET_RULESET 1

/* end private section */





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_UGGRID_CONFIG_HH



#ifndef DUNE_GRID_CONFIG_HH
#define DUNE_GRID_CONFIG_HH

/* Define to 1 if you have module dune-grid available */
#ifndef HAVE_DUNE_GRID
#cmakedefine01 HAVE_DUNE_GRID
#endif




/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "${DUNE_GRID_VERSION}"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR ${DUNE_GRID_VERSION_MAJOR}

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR ${DUNE_GRID_VERSION_MINOR}

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION ${DUNE_GRID_VERSION_REVISION}

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
#cmakedefine DUNE_ALBERTA_VERSION @DUNE_ALBERTA_VERSION@

/* Define to 1 if you have mkstemp function */
#cmakedefine01 HAVE_MKSTEMP







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
#cmakedefine01 HAVE_DUNE_ISTL
#endif





/* Define to the integer type that SuperLU was compiled for
   See e.g. what int_t is defined to in slu_sdefs.h */
#cmakedefine SUPERLU_INT_TYPE @SUPERLU_INT_TYPE@

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "${DUNE_ISTL_VERSION}"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR ${DUNE_ISTL_VERSION_MAJOR}

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR ${DUNE_ISTL_VERSION_MINOR}

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION ${DUNE_ISTL_VERSION_REVISION}

/* Enable/Disable the backwards compatibility of the category enum/method in dune-istl solvers, preconditioner, etc. */
#cmakedefine DUNE_ISTL_SUPPORT_OLD_CATEGORY_INTERFACE @DUNE_ISTL_SUPPORT_OLD_CATEGORY_INTERFACE@





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_ISTL_CONFIG_HH



#ifndef DUNE_TYPETREE_CONFIG_HH
#define DUNE_TYPETREE_CONFIG_HH

/* Define to 1 if you have module dune-typetree available */
#ifndef HAVE_DUNE_TYPETREE
#cmakedefine01 HAVE_DUNE_TYPETREE
#endif




/* Define to the version of dune-typetree */
#define DUNE_TYPETREE_VERSION "${DUNE_TYPETREE_VERSION}"

/* Define to the major version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MAJOR ${DUNE_TYPETREE_VERSION_MAJOR}

/* Define to the minor version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MINOR ${DUNE_TYPETREE_VERSION_MINOR}

/* Define to the revision of dune-typetree */
#define DUNE_TYPETREE_VERSION_REVISION ${DUNE_TYPETREE_VERSION_REVISION}





#if __has_include(<dune-common-config.hh>)
  #include <dune-common-config.hh>
#endif


#endif // DUNE_TYPETREE_CONFIG_HH



#ifndef DUNE_FUNCTIONS_CONFIG_HH
#define DUNE_FUNCTIONS_CONFIG_HH

/* Define to 1 if you have module dune-functions available */
#ifndef HAVE_DUNE_FUNCTIONS
#cmakedefine01 HAVE_DUNE_FUNCTIONS
#endif





/* Define to the version of dune-functions */
#define DUNE_FUNCTIONS_VERSION "@DUNE_FUNCTIONS_VERSION@"

/* Define to the major version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MAJOR @DUNE_FUNCTIONS_VERSION_MAJOR@

/* Define to the minor version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MINOR @DUNE_FUNCTIONS_VERSION_MINOR@

/* Define to the revision of dune-functions */
#define DUNE_FUNCTIONS_VERSION_REVISION @DUNE_FUNCTIONS_VERSION_REVISION@





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
#cmakedefine01 HAVE_DUNE_FOAMGRID
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
#cmakedefine01 HAVE_DUNE_ALUGRID
#endif





#define DUNE_ALUGRID_VERSION "${DUNE_ALUGRID_VERSION}"

/* Define to the major version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MAJOR ${DUNE_ALUGRID_VERSION_MAJOR}

/* Define to the minor version of dune-alugrid */
#define DUNE_ALUGRID_VERSION_MINOR ${DUNE_ALUGRID_VERSION_MINOR}

/* Define to the revision of dune-alugrid*/
#define DUNE_ALUGRID_VERSION_REVISION ${DUNE_ALUGRID_VERSION_REVISION}

/* Define to build more .cc into library */
#cmakedefine DUNE_ALUGRID_COMPILE_BINDINGS_IN_LIB 1

/* Define if we have dlmalloc */
#cmakedefine HAVE_DLMALLOC 1

/* Define if we have zoltan */
#cmakedefine HAVE_ZOLTAN 1

/* Define if we have ZLIB */
#cmakedefine HAVE_ZLIB 1

/* Include source file for dlmalloc */
#cmakedefine DLMALLOC_SOURCE_INCLUDE ${DLMALLOC_SOURCE_INCLUDE}

/* Define if we have thread local storage */
#cmakedefine HAVE_PTHREAD_TLS 1

/* Define if we have pthreads */
#cmakedefine HAVE_PTHREAD 1

/* Define if testgrids.hh from dune-grid have been found in docs/grids/gridfactory */
#cmakedefine HAVE_DUNE_GRID_TESTGRIDS 1

/* Define if METIS is enabled for ALUGrid */
#cmakedefine01 ALUGRID_HAVE_METIS

/* Grid type magic for DGF parser */
@ALUGRID_CONFIG_H_BOTTOM@




#if __has_include(<dune-grid-config.hh>)
  #include <dune-grid-config.hh>
#endif


#endif // DUNE_ALUGRID_CONFIG_HH



#ifndef DUNE_VTK_CONFIG_HH
#define DUNE_VTK_CONFIG_HH

/* Define to 1 if you have module dune-vtk available */
#ifndef HAVE_DUNE_VTK
#cmakedefine01 HAVE_DUNE_VTK
#endif





/* Define to the version of dune-vtk */
#define DUNE_VTK_VERSION "@DUNE_VTK_VERSION@"

/* Define to the major version of dune-vtk */
#define DUNE_VTK_VERSION_MAJOR @DUNE_VTK_VERSION_MAJOR@

/* Define to the minor version of dune-vtk */
#define DUNE_VTK_VERSION_MINOR @DUNE_VTK_VERSION_MINOR@

/* Define to the revision of dune-vtk */
#define DUNE_VTK_VERSION_REVISION @DUNE_VTK_VERSION_REVISION@

/* Define if you have the ZLIB library.  */
#cmakedefine HAVE_VTK_ZLIB ENABLE_VTK_ZLIB





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
#cmakedefine01 HAVE_DUNE_CURVEDGEOMETRY
#endif





/* Define to the version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION "@DUNE_CURVEDGEOMETRY_VERSION@"

/* Define to the major version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_MAJOR @DUNE_CURVEDGEOMETRY_VERSION_MAJOR@

/* Define to the minor version of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_MINOR @DUNE_CURVEDGEOMETRY_VERSION_MINOR@

/* Define to the revision of dune-curvedgeometry */
#define DUNE_CURVEDGEOMETRY_VERSION_REVISION @DUNE_CURVEDGEOMETRY_VERSION_REVISION@





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
#cmakedefine01 HAVE_DUNE_CURVEDGRID
#endif





/* Define to the version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION "@DUNE_CURVEDGRID_VERSION@"

/* Define to the major version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_MAJOR @DUNE_CURVEDGRID_VERSION_MAJOR@

/* Define to the minor version of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_MINOR @DUNE_CURVEDGRID_VERSION_MINOR@

/* Define to the revision of dune-curvedgrid */
#define DUNE_CURVEDGRID_VERSION_REVISION @DUNE_CURVEDGRID_VERSION_REVISION@





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
#cmakedefine01 HAVE_DUNE_MESHDIST
#endif

/* begin hausdorff-distance
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/



/* Define to the version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION "@DUNE_MESHDIST_VERSION@"

/* Define to the major version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_MAJOR @DUNE_MESHDIST_VERSION_MAJOR@

/* Define to the minor version of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_MINOR @DUNE_MESHDIST_VERSION_MINOR@

/* Define to the revision of hausdorff-distance */
#define DUNE_MESHDIST_VERSION_REVISION @DUNE_MESHDIST_VERSION_REVISION@

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
#cmakedefine01 HAVE_DUNE_GMSH4
#endif





/* Define to the version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION "@DUNE_GMSH4_VERSION@"

/* Define to the major version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_MAJOR @DUNE_GMSH4_VERSION_MAJOR@

/* Define to the minor version of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_MINOR @DUNE_GMSH4_VERSION_MINOR@

/* Define to the revision of dune-gmsh4 */
#define DUNE_GMSH4_VERSION_REVISION @DUNE_GMSH4_VERSION_REVISION@





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
#cmakedefine01 HAVE_HIGHER_ORDER_BGN
#endif





/* Define to the version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION "@HIGHER_ORDER_BGN_VERSION@"

/* Define to the major version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_MAJOR @HIGHER_ORDER_BGN_VERSION_MAJOR@

/* Define to the minor version of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_MINOR @HIGHER_ORDER_BGN_VERSION_MINOR@

/* Define to the revision of higher-order-bgn */
#define HIGHER_ORDER_BGN_VERSION_REVISION @HIGHER_ORDER_BGN_VERSION_REVISION@





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
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"



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
@GRID_CONFIG_H_BOTTOM@



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

