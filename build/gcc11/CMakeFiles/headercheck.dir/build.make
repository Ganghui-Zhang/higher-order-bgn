# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/leeloos/dune-test/higher-order-bgn

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/leeloos/dune-test/higher-order-bgn/build/gcc11

# Utility rule file for headercheck.

# Include any custom commands dependencies for this target.
include CMakeFiles/headercheck.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/headercheck.dir/progress.make

CMakeFiles/headercheck:
	/usr/bin/cmake -DENABLE_HEADERCHECK= -P /home/leeloos/dune-test/dune-common/cmake/scripts/FinalizeHeadercheck.cmake

headercheck: CMakeFiles/headercheck
headercheck: CMakeFiles/headercheck.dir/build.make
.PHONY : headercheck

# Rule to build all files generated by this target.
CMakeFiles/headercheck.dir/build: headercheck
.PHONY : CMakeFiles/headercheck.dir/build

CMakeFiles/headercheck.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/headercheck.dir/cmake_clean.cmake
.PHONY : CMakeFiles/headercheck.dir/clean

CMakeFiles/headercheck.dir/depend:
	cd /home/leeloos/dune-test/higher-order-bgn/build/gcc11 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/leeloos/dune-test/higher-order-bgn /home/leeloos/dune-test/higher-order-bgn /home/leeloos/dune-test/higher-order-bgn/build/gcc11 /home/leeloos/dune-test/higher-order-bgn/build/gcc11 /home/leeloos/dune-test/higher-order-bgn/build/gcc11/CMakeFiles/headercheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/headercheck.dir/depend

