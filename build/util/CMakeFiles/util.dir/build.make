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
include util/CMakeFiles/util.dir/depend.make

# Include the progress variables for this target.
include util/CMakeFiles/util.dir/progress.make

# Include the compile flags for this target's objects.
include util/CMakeFiles/util.dir/flags.make

util/CMakeFiles/util.dir/arguments.cc.o: util/CMakeFiles/util.dir/flags.make
util/CMakeFiles/util.dir/arguments.cc.o: ../util/arguments.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object util/CMakeFiles/util.dir/arguments.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/arguments.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/util/arguments.cc

util/CMakeFiles/util.dir/arguments.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/arguments.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/util/arguments.cc > CMakeFiles/util.dir/arguments.cc.i

util/CMakeFiles/util.dir/arguments.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/arguments.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/util/arguments.cc -o CMakeFiles/util.dir/arguments.cc.s

util/CMakeFiles/util.dir/arguments.cc.o.requires:

.PHONY : util/CMakeFiles/util.dir/arguments.cc.o.requires

util/CMakeFiles/util.dir/arguments.cc.o.provides: util/CMakeFiles/util.dir/arguments.cc.o.requires
	$(MAKE) -f util/CMakeFiles/util.dir/build.make util/CMakeFiles/util.dir/arguments.cc.o.provides.build
.PHONY : util/CMakeFiles/util.dir/arguments.cc.o.provides

util/CMakeFiles/util.dir/arguments.cc.o.provides.build: util/CMakeFiles/util.dir/arguments.cc.o


util/CMakeFiles/util.dir/file_system.cc.o: util/CMakeFiles/util.dir/flags.make
util/CMakeFiles/util.dir/file_system.cc.o: ../util/file_system.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object util/CMakeFiles/util.dir/file_system.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/file_system.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/util/file_system.cc

util/CMakeFiles/util.dir/file_system.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/file_system.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/util/file_system.cc > CMakeFiles/util.dir/file_system.cc.i

util/CMakeFiles/util.dir/file_system.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/file_system.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/util/file_system.cc -o CMakeFiles/util.dir/file_system.cc.s

util/CMakeFiles/util.dir/file_system.cc.o.requires:

.PHONY : util/CMakeFiles/util.dir/file_system.cc.o.requires

util/CMakeFiles/util.dir/file_system.cc.o.provides: util/CMakeFiles/util.dir/file_system.cc.o.requires
	$(MAKE) -f util/CMakeFiles/util.dir/build.make util/CMakeFiles/util.dir/file_system.cc.o.provides.build
.PHONY : util/CMakeFiles/util.dir/file_system.cc.o.provides

util/CMakeFiles/util.dir/file_system.cc.o.provides.build: util/CMakeFiles/util.dir/file_system.cc.o


util/CMakeFiles/util.dir/ini_parser.cc.o: util/CMakeFiles/util.dir/flags.make
util/CMakeFiles/util.dir/ini_parser.cc.o: ../util/ini_parser.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object util/CMakeFiles/util.dir/ini_parser.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/ini_parser.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/util/ini_parser.cc

util/CMakeFiles/util.dir/ini_parser.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/ini_parser.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/util/ini_parser.cc > CMakeFiles/util.dir/ini_parser.cc.i

util/CMakeFiles/util.dir/ini_parser.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/ini_parser.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/util/ini_parser.cc -o CMakeFiles/util.dir/ini_parser.cc.s

util/CMakeFiles/util.dir/ini_parser.cc.o.requires:

.PHONY : util/CMakeFiles/util.dir/ini_parser.cc.o.requires

util/CMakeFiles/util.dir/ini_parser.cc.o.provides: util/CMakeFiles/util.dir/ini_parser.cc.o.requires
	$(MAKE) -f util/CMakeFiles/util.dir/build.make util/CMakeFiles/util.dir/ini_parser.cc.o.provides.build
.PHONY : util/CMakeFiles/util.dir/ini_parser.cc.o.provides

util/CMakeFiles/util.dir/ini_parser.cc.o.provides.build: util/CMakeFiles/util.dir/ini_parser.cc.o


util/CMakeFiles/util.dir/system.cc.o: util/CMakeFiles/util.dir/flags.make
util/CMakeFiles/util.dir/system.cc.o: ../util/system.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object util/CMakeFiles/util.dir/system.cc.o"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/util.dir/system.cc.o -c /home/demo/桌面/ImageBasedModellingEdu/util/system.cc

util/CMakeFiles/util.dir/system.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/util.dir/system.cc.i"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/demo/桌面/ImageBasedModellingEdu/util/system.cc > CMakeFiles/util.dir/system.cc.i

util/CMakeFiles/util.dir/system.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/util.dir/system.cc.s"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/demo/桌面/ImageBasedModellingEdu/util/system.cc -o CMakeFiles/util.dir/system.cc.s

util/CMakeFiles/util.dir/system.cc.o.requires:

.PHONY : util/CMakeFiles/util.dir/system.cc.o.requires

util/CMakeFiles/util.dir/system.cc.o.provides: util/CMakeFiles/util.dir/system.cc.o.requires
	$(MAKE) -f util/CMakeFiles/util.dir/build.make util/CMakeFiles/util.dir/system.cc.o.provides.build
.PHONY : util/CMakeFiles/util.dir/system.cc.o.provides

util/CMakeFiles/util.dir/system.cc.o.provides.build: util/CMakeFiles/util.dir/system.cc.o


# Object files for target util
util_OBJECTS = \
"CMakeFiles/util.dir/arguments.cc.o" \
"CMakeFiles/util.dir/file_system.cc.o" \
"CMakeFiles/util.dir/ini_parser.cc.o" \
"CMakeFiles/util.dir/system.cc.o"

# External object files for target util
util_EXTERNAL_OBJECTS =

util/libutil.a: util/CMakeFiles/util.dir/arguments.cc.o
util/libutil.a: util/CMakeFiles/util.dir/file_system.cc.o
util/libutil.a: util/CMakeFiles/util.dir/ini_parser.cc.o
util/libutil.a: util/CMakeFiles/util.dir/system.cc.o
util/libutil.a: util/CMakeFiles/util.dir/build.make
util/libutil.a: util/CMakeFiles/util.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/demo/桌面/ImageBasedModellingEdu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libutil.a"
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && $(CMAKE_COMMAND) -P CMakeFiles/util.dir/cmake_clean_target.cmake
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/util.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
util/CMakeFiles/util.dir/build: util/libutil.a

.PHONY : util/CMakeFiles/util.dir/build

util/CMakeFiles/util.dir/requires: util/CMakeFiles/util.dir/arguments.cc.o.requires
util/CMakeFiles/util.dir/requires: util/CMakeFiles/util.dir/file_system.cc.o.requires
util/CMakeFiles/util.dir/requires: util/CMakeFiles/util.dir/ini_parser.cc.o.requires
util/CMakeFiles/util.dir/requires: util/CMakeFiles/util.dir/system.cc.o.requires

.PHONY : util/CMakeFiles/util.dir/requires

util/CMakeFiles/util.dir/clean:
	cd /home/demo/桌面/ImageBasedModellingEdu/build/util && $(CMAKE_COMMAND) -P CMakeFiles/util.dir/cmake_clean.cmake
.PHONY : util/CMakeFiles/util.dir/clean

util/CMakeFiles/util.dir/depend:
	cd /home/demo/桌面/ImageBasedModellingEdu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/demo/桌面/ImageBasedModellingEdu /home/demo/桌面/ImageBasedModellingEdu/util /home/demo/桌面/ImageBasedModellingEdu/build /home/demo/桌面/ImageBasedModellingEdu/build/util /home/demo/桌面/ImageBasedModellingEdu/build/util/CMakeFiles/util.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : util/CMakeFiles/util.dir/depend

