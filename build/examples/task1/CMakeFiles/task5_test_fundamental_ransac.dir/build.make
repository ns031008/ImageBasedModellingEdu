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
include examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/depend.make

# Include the progress variables for this target.
include examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/progress.make

# Include the compile flags for this target's objects.
include examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/flags.make

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/flags.make
examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o: ../examples/task1/task1-4_test_fundamental_ransac.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/examples/task1/task1-4_test_fundamental_ransac.cc

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/examples/task1/task1-4_test_fundamental_ransac.cc > CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.i

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/examples/task1/task1-4_test_fundamental_ransac.cc -o CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.s

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.requires:

.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.requires

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.provides: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.requires
	$(MAKE) -f examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/build.make examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.provides.build
.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.provides

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.provides.build: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o


# Object files for target task5_test_fundamental_ransac
task5_test_fundamental_ransac_OBJECTS = \
"CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o"

# External object files for target task5_test_fundamental_ransac
task5_test_fundamental_ransac_EXTERNAL_OBJECTS =

examples/task1/task5_test_fundamental_ransac: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o
examples/task1/task5_test_fundamental_ransac: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/build.make
examples/task1/task5_test_fundamental_ransac: sfm/libsfm.a
examples/task1/task5_test_fundamental_ransac: util/libutil.a
examples/task1/task5_test_fundamental_ransac: core/libcore.a
examples/task1/task5_test_fundamental_ransac: features/libfeatures.a
examples/task1/task5_test_fundamental_ransac: core/libcore.a
examples/task1/task5_test_fundamental_ransac: util/libutil.a
examples/task1/task5_test_fundamental_ransac: /usr/lib/x86_64-linux-gnu/libpng.so
examples/task1/task5_test_fundamental_ransac: /usr/lib/x86_64-linux-gnu/libz.so
examples/task1/task5_test_fundamental_ransac: /usr/lib/x86_64-linux-gnu/libjpeg.so
examples/task1/task5_test_fundamental_ransac: /usr/lib/x86_64-linux-gnu/libtiff.so
examples/task1/task5_test_fundamental_ransac: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable task5_test_fundamental_ransac"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/task5_test_fundamental_ransac.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/build: examples/task1/task5_test_fundamental_ransac

.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/build

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/requires: examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/task1-4_test_fundamental_ransac.cc.o.requires

.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/requires

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/clean:
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 && $(CMAKE_COMMAND) -P CMakeFiles/task5_test_fundamental_ransac.dir/cmake_clean.cmake
.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/clean

examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/depend:
	cd /home/demo/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/demo/桌面/ImageBasedModellingEdu /home/demo/桌面/ImageBasedModellingEdu/examples/task1 /home/demo/桌面/ImageBasedModellingEdu/build /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1 /home/demo/桌面/ImageBasedModellingEdu/build/examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/task1/CMakeFiles/task5_test_fundamental_ransac.dir/depend

