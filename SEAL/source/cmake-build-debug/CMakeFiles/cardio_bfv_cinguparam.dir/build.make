# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/cardio_bfv_cinguparam.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cardio_bfv_cinguparam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cardio_bfv_cinguparam.dir/flags.make

CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o: CMakeFiles/cardio_bfv_cinguparam.dir/flags.make
CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o: ../cardio-bfv-opt/cardio.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o -c /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cardio-bfv-opt/cardio.cpp

CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cardio-bfv-opt/cardio.cpp > CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.i

CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cardio-bfv-opt/cardio.cpp -o CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.s

# Object files for target cardio_bfv_cinguparam
cardio_bfv_cinguparam_OBJECTS = \
"CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o"

# External object files for target cardio_bfv_cinguparam
cardio_bfv_cinguparam_EXTERNAL_OBJECTS =

cardio_bfv_cinguparam: CMakeFiles/cardio_bfv_cinguparam.dir/cardio-bfv-opt/cardio.cpp.o
cardio_bfv_cinguparam: CMakeFiles/cardio_bfv_cinguparam.dir/build.make
cardio_bfv_cinguparam: /usr/local/lib/libseal-3.6.a
cardio_bfv_cinguparam: CMakeFiles/cardio_bfv_cinguparam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cardio_bfv_cinguparam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cardio_bfv_cinguparam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cardio_bfv_cinguparam.dir/build: cardio_bfv_cinguparam

.PHONY : CMakeFiles/cardio_bfv_cinguparam.dir/build

CMakeFiles/cardio_bfv_cinguparam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cardio_bfv_cinguparam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cardio_bfv_cinguparam.dir/clean

CMakeFiles/cardio_bfv_cinguparam.dir/depend:
	cd /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles/cardio_bfv_cinguparam.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cardio_bfv_cinguparam.dir/depend

