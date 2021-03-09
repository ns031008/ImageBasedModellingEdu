# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/demo/桌面/ImageBasedModellingEdu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/demo/桌面/ImageBasedModellingEdu/build

# Include any dependencies generated for this target.
include surface/CMakeFiles/surface.dir/depend.make

# Include the progress variables for this target.
include surface/CMakeFiles/surface.dir/progress.make

# Include the compile flags for this target's objects.
include surface/CMakeFiles/surface.dir/flags.make

surface/CMakeFiles/surface.dir/basis_function.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/basis_function.cc.o: ../surface/basis_function.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object surface/CMakeFiles/surface.dir/basis_function.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/basis_function.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/basis_function.cc

surface/CMakeFiles/surface.dir/basis_function.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/basis_function.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/basis_function.cc > CMakeFiles/surface.dir/basis_function.cc.i

surface/CMakeFiles/surface.dir/basis_function.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/basis_function.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/basis_function.cc -o CMakeFiles/surface.dir/basis_function.cc.s

surface/CMakeFiles/surface.dir/basis_function.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/basis_function.cc.o.requires

surface/CMakeFiles/surface.dir/basis_function.cc.o.provides: surface/CMakeFiles/surface.dir/basis_function.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/basis_function.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/basis_function.cc.o.provides

surface/CMakeFiles/surface.dir/basis_function.cc.o.provides.build: surface/CMakeFiles/surface.dir/basis_function.cc.o


surface/CMakeFiles/surface.dir/hermite.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/hermite.cc.o: ../surface/hermite.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object surface/CMakeFiles/surface.dir/hermite.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/hermite.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/hermite.cc

surface/CMakeFiles/surface.dir/hermite.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/hermite.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/hermite.cc > CMakeFiles/surface.dir/hermite.cc.i

surface/CMakeFiles/surface.dir/hermite.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/hermite.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/hermite.cc -o CMakeFiles/surface.dir/hermite.cc.s

surface/CMakeFiles/surface.dir/hermite.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/hermite.cc.o.requires

surface/CMakeFiles/surface.dir/hermite.cc.o.provides: surface/CMakeFiles/surface.dir/hermite.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/hermite.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/hermite.cc.o.provides

surface/CMakeFiles/surface.dir/hermite.cc.o.provides.build: surface/CMakeFiles/surface.dir/hermite.cc.o


surface/CMakeFiles/surface.dir/iso_octree.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/iso_octree.cc.o: ../surface/iso_octree.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object surface/CMakeFiles/surface.dir/iso_octree.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/iso_octree.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/iso_octree.cc

surface/CMakeFiles/surface.dir/iso_octree.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/iso_octree.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/iso_octree.cc > CMakeFiles/surface.dir/iso_octree.cc.i

surface/CMakeFiles/surface.dir/iso_octree.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/iso_octree.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/iso_octree.cc -o CMakeFiles/surface.dir/iso_octree.cc.s

surface/CMakeFiles/surface.dir/iso_octree.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/iso_octree.cc.o.requires

surface/CMakeFiles/surface.dir/iso_octree.cc.o.provides: surface/CMakeFiles/surface.dir/iso_octree.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/iso_octree.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/iso_octree.cc.o.provides

surface/CMakeFiles/surface.dir/iso_octree.cc.o.provides.build: surface/CMakeFiles/surface.dir/iso_octree.cc.o


surface/CMakeFiles/surface.dir/iso_surface.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/iso_surface.cc.o: ../surface/iso_surface.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object surface/CMakeFiles/surface.dir/iso_surface.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/iso_surface.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/iso_surface.cc

surface/CMakeFiles/surface.dir/iso_surface.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/iso_surface.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/iso_surface.cc > CMakeFiles/surface.dir/iso_surface.cc.i

surface/CMakeFiles/surface.dir/iso_surface.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/iso_surface.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/iso_surface.cc -o CMakeFiles/surface.dir/iso_surface.cc.s

surface/CMakeFiles/surface.dir/iso_surface.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/iso_surface.cc.o.requires

surface/CMakeFiles/surface.dir/iso_surface.cc.o.provides: surface/CMakeFiles/surface.dir/iso_surface.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/iso_surface.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/iso_surface.cc.o.provides

surface/CMakeFiles/surface.dir/iso_surface.cc.o.provides.build: surface/CMakeFiles/surface.dir/iso_surface.cc.o


surface/CMakeFiles/surface.dir/mesh_clean.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/mesh_clean.cc.o: ../surface/mesh_clean.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object surface/CMakeFiles/surface.dir/mesh_clean.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/mesh_clean.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/mesh_clean.cc

surface/CMakeFiles/surface.dir/mesh_clean.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/mesh_clean.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/mesh_clean.cc > CMakeFiles/surface.dir/mesh_clean.cc.i

surface/CMakeFiles/surface.dir/mesh_clean.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/mesh_clean.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/mesh_clean.cc -o CMakeFiles/surface.dir/mesh_clean.cc.s

surface/CMakeFiles/surface.dir/mesh_clean.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/mesh_clean.cc.o.requires

surface/CMakeFiles/surface.dir/mesh_clean.cc.o.provides: surface/CMakeFiles/surface.dir/mesh_clean.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/mesh_clean.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/mesh_clean.cc.o.provides

surface/CMakeFiles/surface.dir/mesh_clean.cc.o.provides.build: surface/CMakeFiles/surface.dir/mesh_clean.cc.o


surface/CMakeFiles/surface.dir/octree.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/octree.cc.o: ../surface/octree.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object surface/CMakeFiles/surface.dir/octree.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/octree.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/octree.cc

surface/CMakeFiles/surface.dir/octree.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/octree.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/octree.cc > CMakeFiles/surface.dir/octree.cc.i

surface/CMakeFiles/surface.dir/octree.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/octree.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/octree.cc -o CMakeFiles/surface.dir/octree.cc.s

surface/CMakeFiles/surface.dir/octree.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/octree.cc.o.requires

surface/CMakeFiles/surface.dir/octree.cc.o.provides: surface/CMakeFiles/surface.dir/octree.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/octree.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/octree.cc.o.provides

surface/CMakeFiles/surface.dir/octree.cc.o.provides.build: surface/CMakeFiles/surface.dir/octree.cc.o


surface/CMakeFiles/surface.dir/sample_io.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/sample_io.cc.o: ../surface/sample_io.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object surface/CMakeFiles/surface.dir/sample_io.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/sample_io.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/sample_io.cc

surface/CMakeFiles/surface.dir/sample_io.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/sample_io.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/sample_io.cc > CMakeFiles/surface.dir/sample_io.cc.i

surface/CMakeFiles/surface.dir/sample_io.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/sample_io.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/sample_io.cc -o CMakeFiles/surface.dir/sample_io.cc.s

surface/CMakeFiles/surface.dir/sample_io.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/sample_io.cc.o.requires

surface/CMakeFiles/surface.dir/sample_io.cc.o.provides: surface/CMakeFiles/surface.dir/sample_io.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/sample_io.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/sample_io.cc.o.provides

surface/CMakeFiles/surface.dir/sample_io.cc.o.provides.build: surface/CMakeFiles/surface.dir/sample_io.cc.o


surface/CMakeFiles/surface.dir/triangulation.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/triangulation.cc.o: ../surface/triangulation.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object surface/CMakeFiles/surface.dir/triangulation.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/triangulation.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/triangulation.cc

surface/CMakeFiles/surface.dir/triangulation.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/triangulation.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/triangulation.cc > CMakeFiles/surface.dir/triangulation.cc.i

surface/CMakeFiles/surface.dir/triangulation.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/triangulation.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/triangulation.cc -o CMakeFiles/surface.dir/triangulation.cc.s

surface/CMakeFiles/surface.dir/triangulation.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/triangulation.cc.o.requires

surface/CMakeFiles/surface.dir/triangulation.cc.o.provides: surface/CMakeFiles/surface.dir/triangulation.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/triangulation.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/triangulation.cc.o.provides

surface/CMakeFiles/surface.dir/triangulation.cc.o.provides.build: surface/CMakeFiles/surface.dir/triangulation.cc.o


surface/CMakeFiles/surface.dir/voxel.cc.o: surface/CMakeFiles/surface.dir/flags.make
surface/CMakeFiles/surface.dir/voxel.cc.o: ../surface/voxel.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object surface/CMakeFiles/surface.dir/voxel.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/surface.dir/voxel.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/surface/voxel.cc

surface/CMakeFiles/surface.dir/voxel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/surface.dir/voxel.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/surface/voxel.cc > CMakeFiles/surface.dir/voxel.cc.i

surface/CMakeFiles/surface.dir/voxel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/surface.dir/voxel.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/surface/voxel.cc -o CMakeFiles/surface.dir/voxel.cc.s

surface/CMakeFiles/surface.dir/voxel.cc.o.requires:

.PHONY : surface/CMakeFiles/surface.dir/voxel.cc.o.requires

surface/CMakeFiles/surface.dir/voxel.cc.o.provides: surface/CMakeFiles/surface.dir/voxel.cc.o.requires
	$(MAKE) -f surface/CMakeFiles/surface.dir/build.make surface/CMakeFiles/surface.dir/voxel.cc.o.provides.build
.PHONY : surface/CMakeFiles/surface.dir/voxel.cc.o.provides

surface/CMakeFiles/surface.dir/voxel.cc.o.provides.build: surface/CMakeFiles/surface.dir/voxel.cc.o


# Object files for target surface
surface_OBJECTS = \
"CMakeFiles/surface.dir/basis_function.cc.o" \
"CMakeFiles/surface.dir/hermite.cc.o" \
"CMakeFiles/surface.dir/iso_octree.cc.o" \
"CMakeFiles/surface.dir/iso_surface.cc.o" \
"CMakeFiles/surface.dir/mesh_clean.cc.o" \
"CMakeFiles/surface.dir/octree.cc.o" \
"CMakeFiles/surface.dir/sample_io.cc.o" \
"CMakeFiles/surface.dir/triangulation.cc.o" \
"CMakeFiles/surface.dir/voxel.cc.o"

# External object files for target surface
surface_EXTERNAL_OBJECTS =

surface/libsurface.a: surface/CMakeFiles/surface.dir/basis_function.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/hermite.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/iso_octree.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/iso_surface.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/mesh_clean.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/octree.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/sample_io.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/triangulation.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/voxel.cc.o
surface/libsurface.a: surface/CMakeFiles/surface.dir/build.make
surface/libsurface.a: surface/CMakeFiles/surface.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX static library libsurface.a"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && $(CMAKE_COMMAND) -P CMakeFiles/surface.dir/cmake_clean_target.cmake
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/surface.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
surface/CMakeFiles/surface.dir/build: surface/libsurface.a

.PHONY : surface/CMakeFiles/surface.dir/build

surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/basis_function.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/hermite.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/iso_octree.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/iso_surface.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/mesh_clean.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/octree.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/sample_io.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/triangulation.cc.o.requires
surface/CMakeFiles/surface.dir/requires: surface/CMakeFiles/surface.dir/voxel.cc.o.requires

.PHONY : surface/CMakeFiles/surface.dir/requires

surface/CMakeFiles/surface.dir/clean:
	cd /home/demo/桌面/ImageBasedModellingEdu/build/surface && $(CMAKE_COMMAND) -P CMakeFiles/surface.dir/cmake_clean.cmake
.PHONY : surface/CMakeFiles/surface.dir/clean

surface/CMakeFiles/surface.dir/depend:
	cd /home/demo/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/demo/桌面/ImageBasedModellingEdu /home/demo/桌面/ImageBasedModellingEdu/surface /home/demo/桌面/ImageBasedModellingEdu/build /home/demo/桌面/ImageBasedModellingEdu/build/surface /home/demo/桌面/ImageBasedModellingEdu/build/surface/CMakeFiles/surface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : surface/CMakeFiles/surface.dir/depend

