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
include examples/task2/CMakeFiles/task2-1_test_triangle.dir/depend.make

# Include the progress variables for this target.
include examples/task2/CMakeFiles/task2-1_test_triangle.dir/progress.make

# Include the compile flags for this target's objects.
include examples/task2/CMakeFiles/task2-1_test_triangle.dir/flags.make

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o: examples/task2/CMakeFiles/task2-1_test_triangle.dir/flags.make
examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o: ../examples/task2/task2-1_test_triangle.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/examples/task2/task2-1_test_triangle.cc

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/examples/task2/task2-1_test_triangle.cc > CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.i

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/examples/task2/task2-1_test_triangle.cc -o CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.s

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.requires:

.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.requires

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.provides: examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.requires
	$(MAKE) -f examples/task2/CMakeFiles/task2-1_test_triangle.dir/build.make examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.provides.build
.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.provides

examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.provides.build: examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o


# Object files for target task2-1_test_triangle
task2__1_test_triangle_OBJECTS = \
"CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o"

# External object files for target task2-1_test_triangle
task2__1_test_triangle_EXTERNAL_OBJECTS =

examples/task2/task2-1_test_triangle: examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o
examples/task2/task2-1_test_triangle: examples/task2/CMakeFiles/task2-1_test_triangle.dir/build.make
examples/task2/task2-1_test_triangle: sfm/libsfm.a
examples/task2/task2-1_test_triangle: util/libutil.a
examples/task2/task2-1_test_triangle: core/libcore.a
examples/task2/task2-1_test_triangle: util/libutil.a
examples/task2/task2-1_test_triangle: /usr/lib/x86_64-linux-gnu/libpng.so
examples/task2/task2-1_test_triangle: /usr/lib/x86_64-linux-gnu/libz.so
examples/task2/task2-1_test_triangle: /usr/lib/x86_64-linux-gnu/libjpeg.so
examples/task2/task2-1_test_triangle: /usr/lib/x86_64-linux-gnu/libtiff.so
examples/task2/task2-1_test_triangle: examples/task2/CMakeFiles/task2-1_test_triangle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable task2-1_test_triangle"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/task2-1_test_triangle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/task2/CMakeFiles/task2-1_test_triangle.dir/build: examples/task2/task2-1_test_triangle

.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/build

examples/task2/CMakeFiles/task2-1_test_triangle.dir/requires: examples/task2/CMakeFiles/task2-1_test_triangle.dir/task2-1_test_triangle.cc.o.requires

.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/requires

examples/task2/CMakeFiles/task2-1_test_triangle.dir/clean:
	cd /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 && $(CMAKE_COMMAND) -P CMakeFiles/task2-1_test_triangle.dir/cmake_clean.cmake
.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/clean

examples/task2/CMakeFiles/task2-1_test_triangle.dir/depend:
	cd /home/demo/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/demo/桌面/ImageBasedModellingEdu /home/demo/桌面/ImageBasedModellingEdu/examples/task2 /home/demo/桌面/ImageBasedModellingEdu/build /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2 /home/demo/桌面/ImageBasedModellingEdu/build/examples/task2/CMakeFiles/task2-1_test_triangle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/task2/CMakeFiles/task2-1_test_triangle.dir/depend

