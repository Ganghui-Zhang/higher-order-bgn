# Install script for directory: /home/leeloos/dune-test/higher-order-bgn/doc/doxygen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND /usr/bin/cmake --build /home/leeloos/dune-test/higher-order-bgn/build/gcc11 --target doxygen_higher-order-bgn
          WORKING_DIRECTORY /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen)
        file(GLOB doxygenfiles
          GLOB /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/html/*.html
          /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/html/*.js
          /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/html/*.png
          /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/html/*.css
          /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/html/*.gif
          /home/leeloos/dune-test/higher-order-bgn/build/gcc11/doc/doxygen/*.tag
          )
        set(doxygenfiles "${doxygenfiles}")
        foreach(_file ${doxygenfiles})
           get_filename_component(_basename ${_file} NAME)
           # Manifest is generated when prefix was set at configuration time, otherwise is skipped
           LIST(APPEND CMAKE_INSTALL_MANIFEST_FILES /usr/local/share/doc/higher-order-bgn/doxygen/${_basename})
         endforeach()
         file(INSTALL ${doxygenfiles} DESTINATION share/doc/higher-order-bgn/doxygen)
         message(STATUS "Installed doxygen into share/doc/higher-order-bgn/doxygen")
endif()

