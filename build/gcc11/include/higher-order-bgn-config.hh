
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


#ifndef HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH
#define HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH



#endif // HIGHER_ORDER_BGN_CONFIG_BOTTOM_HH
