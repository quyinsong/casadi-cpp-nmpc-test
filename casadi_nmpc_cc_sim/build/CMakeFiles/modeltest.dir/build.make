# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/qys/mySim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qys/mySim/build

# Include any dependencies generated for this target.
include CMakeFiles/modeltest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/modeltest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/modeltest.dir/flags.make

CMakeFiles/modeltest.dir/test/main1.cpp.o: CMakeFiles/modeltest.dir/flags.make
CMakeFiles/modeltest.dir/test/main1.cpp.o: ../test/main1.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qys/mySim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/modeltest.dir/test/main1.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/modeltest.dir/test/main1.cpp.o -c /home/qys/mySim/test/main1.cpp

CMakeFiles/modeltest.dir/test/main1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/modeltest.dir/test/main1.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qys/mySim/test/main1.cpp > CMakeFiles/modeltest.dir/test/main1.cpp.i

CMakeFiles/modeltest.dir/test/main1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/modeltest.dir/test/main1.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qys/mySim/test/main1.cpp -o CMakeFiles/modeltest.dir/test/main1.cpp.s

# Object files for target modeltest
modeltest_OBJECTS = \
"CMakeFiles/modeltest.dir/test/main1.cpp.o"

# External object files for target modeltest
modeltest_EXTERNAL_OBJECTS =

modeltest: CMakeFiles/modeltest.dir/test/main1.cpp.o
modeltest: CMakeFiles/modeltest.dir/build.make
modeltest: CMakeFiles/modeltest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qys/mySim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable modeltest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/modeltest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/modeltest.dir/build: modeltest

.PHONY : CMakeFiles/modeltest.dir/build

CMakeFiles/modeltest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/modeltest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/modeltest.dir/clean

CMakeFiles/modeltest.dir/depend:
	cd /home/qys/mySim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qys/mySim /home/qys/mySim /home/qys/mySim/build /home/qys/mySim/build /home/qys/mySim/build/CMakeFiles/modeltest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/modeltest.dir/depend
