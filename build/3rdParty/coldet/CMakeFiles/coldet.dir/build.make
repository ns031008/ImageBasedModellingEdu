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
include 3rdParty/coldet/CMakeFiles/coldet.dir/depend.make

# Include the progress variables for this target.
include 3rdParty/coldet/CMakeFiles/coldet.dir/progress.make

# Include the compile flags for this target's objects.
include 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o: ../3rdParty/coldet/src/box_bld.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/box_bld.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box_bld.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/box_bld.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box_bld.cpp > CMakeFiles/coldet.dir/src/box_bld.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/box_bld.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box_bld.cpp -o CMakeFiles/coldet.dir/src/box_bld.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o: ../3rdParty/coldet/src/box.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/box.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/box.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box.cpp > CMakeFiles/coldet.dir/src/box.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/box.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/box.cpp -o CMakeFiles/coldet.dir/src/box.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o: ../3rdParty/coldet/src/cdmath3d.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/cdmath3d.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/cdmath3d.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/cdmath3d.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/cdmath3d.cpp > CMakeFiles/coldet.dir/src/cdmath3d.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/cdmath3d.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/cdmath3d.cpp -o CMakeFiles/coldet.dir/src/cdmath3d.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o: ../3rdParty/coldet/src/coldet_bld.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/coldet_bld.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet_bld.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/coldet_bld.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet_bld.cpp > CMakeFiles/coldet.dir/src/coldet_bld.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/coldet_bld.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet_bld.cpp -o CMakeFiles/coldet.dir/src/coldet_bld.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o: ../3rdParty/coldet/src/coldet.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/coldet.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/coldet.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet.cpp > CMakeFiles/coldet.dir/src/coldet.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/coldet.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/coldet.cpp -o CMakeFiles/coldet.dir/src/coldet.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o: ../3rdParty/coldet/src/multiobject.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/multiobject.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/multiobject.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/multiobject.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/multiobject.cpp > CMakeFiles/coldet.dir/src/multiobject.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/multiobject.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/multiobject.cpp -o CMakeFiles/coldet.dir/src/multiobject.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o: ../3rdParty/coldet/src/mytritri.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/mytritri.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/mytritri.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/mytritri.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/mytritri.cpp > CMakeFiles/coldet.dir/src/mytritri.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/mytritri.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/mytritri.cpp -o CMakeFiles/coldet.dir/src/mytritri.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o: ../3rdParty/coldet/src/sysdep.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/coldet.dir/src/sysdep.cpp.o -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/sysdep.cpp

3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coldet.dir/src/sysdep.cpp.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/sysdep.cpp > CMakeFiles/coldet.dir/src/sysdep.cpp.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coldet.dir/src/sysdep.cpp.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/sysdep.cpp -o CMakeFiles/coldet.dir/src/sysdep.cpp.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o


3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o: 3rdParty/coldet/CMakeFiles/coldet.dir/flags.make
3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o: ../3rdParty/coldet/src/tritri.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/gcc-5  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/coldet.dir/src/tritri.c.o   -c /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/tritri.c

3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/coldet.dir/src/tritri.c.i"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/gcc-5  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/tritri.c > CMakeFiles/coldet.dir/src/tritri.c.i

3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/coldet.dir/src/tritri.c.s"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && /usr/bin/gcc-5  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet/src/tritri.c -o CMakeFiles/coldet.dir/src/tritri.c.s

3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.requires:

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.requires

3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.provides: 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.requires
	$(MAKE) -f 3rdParty/coldet/CMakeFiles/coldet.dir/build.make 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.provides.build
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.provides

3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.provides.build: 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o


# Object files for target coldet
coldet_OBJECTS = \
"CMakeFiles/coldet.dir/src/box_bld.cpp.o" \
"CMakeFiles/coldet.dir/src/box.cpp.o" \
"CMakeFiles/coldet.dir/src/cdmath3d.cpp.o" \
"CMakeFiles/coldet.dir/src/coldet_bld.cpp.o" \
"CMakeFiles/coldet.dir/src/coldet.cpp.o" \
"CMakeFiles/coldet.dir/src/multiobject.cpp.o" \
"CMakeFiles/coldet.dir/src/mytritri.cpp.o" \
"CMakeFiles/coldet.dir/src/sysdep.cpp.o" \
"CMakeFiles/coldet.dir/src/tritri.c.o"

# External object files for target coldet
coldet_EXTERNAL_OBJECTS =

3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/build.make
3rdParty/coldet/libcoldet.a: 3rdParty/coldet/CMakeFiles/coldet.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chen/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX static library libcoldet.a"
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && $(CMAKE_COMMAND) -P CMakeFiles/coldet.dir/cmake_clean_target.cmake
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/coldet.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
3rdParty/coldet/CMakeFiles/coldet.dir/build: 3rdParty/coldet/libcoldet.a

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/build

3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box_bld.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/box.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/cdmath3d.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet_bld.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/coldet.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/multiobject.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/mytritri.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/sysdep.cpp.o.requires
3rdParty/coldet/CMakeFiles/coldet.dir/requires: 3rdParty/coldet/CMakeFiles/coldet.dir/src/tritri.c.o.requires

.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/requires

3rdParty/coldet/CMakeFiles/coldet.dir/clean:
	cd /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet && $(CMAKE_COMMAND) -P CMakeFiles/coldet.dir/cmake_clean.cmake
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/clean

3rdParty/coldet/CMakeFiles/coldet.dir/depend:
	cd /home/chen/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chen/桌面/ImageBasedModellingEdu /home/chen/桌面/ImageBasedModellingEdu/3rdParty/coldet /home/chen/桌面/ImageBasedModellingEdu/build /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet /home/chen/桌面/ImageBasedModellingEdu/build/3rdParty/coldet/CMakeFiles/coldet.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : 3rdParty/coldet/CMakeFiles/coldet.dir/depend

