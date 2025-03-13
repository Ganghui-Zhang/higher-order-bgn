cmake_minimum_required(VERSION 3.16)
project(higher-order-bgn CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# deactivate global include-directories (old versions of dune used to need the OLD policy to work)
dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)

# start a dune project with information from dune.module
dune_project()

add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
