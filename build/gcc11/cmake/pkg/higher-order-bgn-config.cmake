if(NOT higher-order-bgn_FOUND)
# Whether this module is installed or not
set(higher-order-bgn_INSTALLED ON)

# Settings specific to the module

# Package initialization

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was higher-order-bgn-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#report other information
set_and_check(higher-order-bgn_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(higher-order-bgn_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set(higher-order-bgn_CMAKE_CONFIG_VERSION "2.11-git")
set(higher-order-bgn_CXX_FLAGS "-march=native")
set(higher-order-bgn_CXX_FLAGS_DEBUG "-g")
set(higher-order-bgn_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(higher-order-bgn_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(higher-order-bgn_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(higher-order-bgn_DEPENDS "dune-curvedgrid;dune-meshdist;dune-gmsh4;dune-vtk;dune-alugrid;dune-foamgrid")
set(higher-order-bgn_SUGGESTS "")
set(higher-order-bgn_MODULE_PATH "${PACKAGE_PREFIX_DIR}/share/dune/cmake/modules")
set(higher-order-bgn_PYTHON_WHEELHOUSE "${PACKAGE_PREFIX_DIR}/share/dune/wheelhouse")
set(higher-order-bgn_LIBRARIES "")
set(higher-order-bgn_HASPYTHON 0)
set(higher-order-bgn_PYTHONREQUIRES "")

# Resolve dune dependencies
include(CMakeFindDependencyMacro)
macro(find_and_check_dune_dependency module version)
  find_dependency(${module})
  list(PREPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
  include(DuneModuleDependencies)
  list(POP_FRONT CMAKE_MODULE_PATH)
  if(dune-common_VERSION VERSION_GREATER_EQUAL "2.10")
    dune_check_module_version(${module} QUIET REQUIRED VERSION "${version}")
  endif()
endmacro()

find_and_check_dune_dependency(dune-curvedgrid " ")
find_and_check_dune_dependency(dune-meshdist " ")
find_and_check_dune_dependency(dune-gmsh4 " ")
find_and_check_dune_dependency(dune-vtk " ")
find_and_check_dune_dependency(dune-alugrid " ")
find_and_check_dune_dependency(dune-foamgrid " ")

# Set up DUNE_LIBS, DUNE_FOUND_DEPENDENCIES, DUNE_*_FOUND, and HAVE_* variables
if(higher-order-bgn_LIBRARIES)
  message(STATUS "Setting higher-order-bgn_LIBRARIES=${higher-order-bgn_LIBRARIES}")
  list(PREPEND DUNE_LIBS ${higher-order-bgn_LIBRARIES})
endif()
list(APPEND DUNE_FOUND_DEPENDENCIES higher-order-bgn)
set(DUNE_higher-order-bgn_FOUND TRUE)
set(HAVE_HIGHER_ORDER_BGN TRUE)

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


# If this file is found in a super build that includes higher-order-bgn, the
# `higher-order-bgn-targets.cmake`-file has not yet been generated. This variable
# determines whether the configuration of higher-order-bgn has been completed.
get_property(higher-order-bgn_IN_CONFIG_MODE GLOBAL PROPERTY higher-order-bgn_LIBRARIES DEFINED)

#import the target
if(higher-order-bgn_LIBRARIES AND NOT higher-order-bgn_IN_CONFIG_MODE)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/higher-order-bgn-targets.cmake")
endif()

endif()
