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
CMAKE_SOURCE_DIR = /home/chen/桌面/ImageBasedModellingEdu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chen/桌面/ImageBasedModellingEdu/build

# Include any dependencies generated for this target.
include mvs/CMakeFiles/mvs.dir/depend.make

# Include the progress variables for this target.
include mvs/CMakeFiles/mvs.dir/progress.make

# Include the compile flags for this target's objects.
include mvs/CMakeFiles/mvs.dir/flags.make

mvs/CMakeFiles/mvs.dir/dmrecon.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/dmrecon.cc.o: ../mvs/dmrecon.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object mvs/CMakeFiles/mvs.dir/dmrecon.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/dmrecon.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/dmrecon.cc

mvs/CMakeFiles/mvs.dir/dmrecon.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/dmrecon.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/dmrecon.cc > CMakeFiles/mvs.dir/dmrecon.cc.i

mvs/CMakeFiles/mvs.dir/dmrecon.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/dmrecon.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/dmrecon.cc -o CMakeFiles/mvs.dir/dmrecon.cc.s

mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.requires

mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.provides: mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.provides

mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/dmrecon.cc.o


mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o: ../mvs/global_view_selection.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/global_view_selection.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/global_view_selection.cc

mvs/CMakeFiles/mvs.dir/global_view_selection.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/global_view_selection.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/global_view_selection.cc > CMakeFiles/mvs.dir/global_view_selection.cc.i

mvs/CMakeFiles/mvs.dir/global_view_selection.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/global_view_selection.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/global_view_selection.cc -o CMakeFiles/mvs.dir/global_view_selection.cc.s

mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.requires

mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.provides: mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.provides

mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o


mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o: ../mvs/image_pyramid.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/image_pyramid.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/image_pyramid.cc

mvs/CMakeFiles/mvs.dir/image_pyramid.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/image_pyramid.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/image_pyramid.cc > CMakeFiles/mvs.dir/image_pyramid.cc.i

mvs/CMakeFiles/mvs.dir/image_pyramid.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/image_pyramid.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/image_pyramid.cc -o CMakeFiles/mvs.dir/image_pyramid.cc.s

mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.requires

mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.provides: mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.provides

mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o


mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o: ../mvs/local_view_selection.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/local_view_selection.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/local_view_selection.cc

mvs/CMakeFiles/mvs.dir/local_view_selection.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/local_view_selection.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/local_view_selection.cc > CMakeFiles/mvs.dir/local_view_selection.cc.i

mvs/CMakeFiles/mvs.dir/local_view_selection.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/local_view_selection.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/local_view_selection.cc -o CMakeFiles/mvs.dir/local_view_selection.cc.s

mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.requires

mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.provides: mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.provides

mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o


mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o: ../mvs/mvs_tools.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/mvs_tools.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/mvs_tools.cc

mvs/CMakeFiles/mvs.dir/mvs_tools.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/mvs_tools.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/mvs_tools.cc > CMakeFiles/mvs.dir/mvs_tools.cc.i

mvs/CMakeFiles/mvs.dir/mvs_tools.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/mvs_tools.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/mvs_tools.cc -o CMakeFiles/mvs.dir/mvs_tools.cc.s

mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.requires

mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.provides: mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.provides

mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o


mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o: ../mvs/patch_optimization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/patch_optimization.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_optimization.cc

mvs/CMakeFiles/mvs.dir/patch_optimization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/patch_optimization.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_optimization.cc > CMakeFiles/mvs.dir/patch_optimization.cc.i

mvs/CMakeFiles/mvs.dir/patch_optimization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/patch_optimization.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_optimization.cc -o CMakeFiles/mvs.dir/patch_optimization.cc.s

mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.requires

mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.provides: mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.provides

mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o


mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o: ../mvs/patch_sampler.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/patch_sampler.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_sampler.cc

mvs/CMakeFiles/mvs.dir/patch_sampler.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/patch_sampler.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_sampler.cc > CMakeFiles/mvs.dir/patch_sampler.cc.i

mvs/CMakeFiles/mvs.dir/patch_sampler.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/patch_sampler.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/patch_sampler.cc -o CMakeFiles/mvs.dir/patch_sampler.cc.s

mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.requires

mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.provides: mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.provides

mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o


mvs/CMakeFiles/mvs.dir/single_view.cc.o: mvs/CMakeFiles/mvs.dir/flags.make
mvs/CMakeFiles/mvs.dir/single_view.cc.o: ../mvs/single_view.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object mvs/CMakeFiles/mvs.dir/single_view.cc.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mvs.dir/single_view.cc.o -c /home/chen/桌面/ImageBasedModellingEdu/mvs/single_view.cc

mvs/CMakeFiles/mvs.dir/single_view.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mvs.dir/single_view.cc.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/mvs/single_view.cc > CMakeFiles/mvs.dir/single_view.cc.i

mvs/CMakeFiles/mvs.dir/single_view.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mvs.dir/single_view.cc.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/mvs/single_view.cc -o CMakeFiles/mvs.dir/single_view.cc.s

mvs/CMakeFiles/mvs.dir/single_view.cc.o.requires:

.PHONY : mvs/CMakeFiles/mvs.dir/single_view.cc.o.requires

mvs/CMakeFiles/mvs.dir/single_view.cc.o.provides: mvs/CMakeFiles/mvs.dir/single_view.cc.o.requires
	$(MAKE) -f mvs/CMakeFiles/mvs.dir/build.make mvs/CMakeFiles/mvs.dir/single_view.cc.o.provides.build
.PHONY : mvs/CMakeFiles/mvs.dir/single_view.cc.o.provides

mvs/CMakeFiles/mvs.dir/single_view.cc.o.provides.build: mvs/CMakeFiles/mvs.dir/single_view.cc.o


# Object files for target mvs
mvs_OBJECTS = \
"CMakeFiles/mvs.dir/dmrecon.cc.o" \
"CMakeFiles/mvs.dir/global_view_selection.cc.o" \
"CMakeFiles/mvs.dir/image_pyramid.cc.o" \
"CMakeFiles/mvs.dir/local_view_selection.cc.o" \
"CMakeFiles/mvs.dir/mvs_tools.cc.o" \
"CMakeFiles/mvs.dir/patch_optimization.cc.o" \
"CMakeFiles/mvs.dir/patch_sampler.cc.o" \
"CMakeFiles/mvs.dir/single_view.cc.o"

# External object files for target mvs
mvs_EXTERNAL_OBJECTS =

mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/dmrecon.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/single_view.cc.o
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/build.make
mvs/libmvs.a: mvs/CMakeFiles/mvs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX static library libmvs.a"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && $(CMAKE_COMMAND) -P CMakeFiles/mvs.dir/cmake_clean_target.cmake
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mvs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
mvs/CMakeFiles/mvs.dir/build: mvs/libmvs.a

.PHONY : mvs/CMakeFiles/mvs.dir/build

mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/dmrecon.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/global_view_selection.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/image_pyramid.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/local_view_selection.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/mvs_tools.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/patch_optimization.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/patch_sampler.cc.o.requires
mvs/CMakeFiles/mvs.dir/requires: mvs/CMakeFiles/mvs.dir/single_view.cc.o.requires

.PHONY : mvs/CMakeFiles/mvs.dir/requires

mvs/CMakeFiles/mvs.dir/clean:
	cd /home/chen/桌面/ImageBasedModellingEdu/build/mvs && $(CMAKE_COMMAND) -P CMakeFiles/mvs.dir/cmake_clean.cmake
.PHONY : mvs/CMakeFiles/mvs.dir/clean

mvs/CMakeFiles/mvs.dir/depend:
	cd /home/chen/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chen/桌面/ImageBasedModellingEdu /home/chen/桌面/ImageBasedModellingEdu/mvs /home/chen/桌面/ImageBasedModellingEdu/build /home/chen/桌面/ImageBasedModellingEdu/build/mvs /home/chen/桌面/ImageBasedModellingEdu/build/mvs/CMakeFiles/mvs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : mvs/CMakeFiles/mvs.dir/depend

