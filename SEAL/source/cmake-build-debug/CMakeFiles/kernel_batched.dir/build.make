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
include CMakeFiles/kernel_batched.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/kernel_batched.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/kernel_batched.dir/flags.make

CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o: CMakeFiles/kernel_batched.dir/flags.make
CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o: ../kernel-bfv-batched/kernel_batched.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o -c /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/kernel-bfv-batched/kernel_batched.cpp

CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/kernel-bfv-batched/kernel_batched.cpp > CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.i

CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/kernel-bfv-batched/kernel_batched.cpp -o CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.s

# Object files for target kernel_batched
kernel_batched_OBJECTS = \
"CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o"

# External object files for target kernel_batched
kernel_batched_EXTERNAL_OBJECTS =

kernel_batched: CMakeFiles/kernel_batched.dir/kernel-bfv-batched/kernel_batched.cpp.o
kernel_batched: CMakeFiles/kernel_batched.dir/build.make
kernel_batched: /usr/local/lib/libseal-3.6.a
kernel_batched: CMakeFiles/kernel_batched.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable kernel_batched"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kernel_batched.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/kernel_batched.dir/build: kernel_batched

.PHONY : CMakeFiles/kernel_batched.dir/build

CMakeFiles/kernel_batched.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/kernel_batched.dir/cmake_clean.cmake
.PHONY : CMakeFiles/kernel_batched.dir/clean

CMakeFiles/kernel_batched.dir/depend:
	cd /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug /Users/mwinger/Desktop/MasterThesis/SoK/SEAL/source/cmake-build-debug/CMakeFiles/kernel_batched.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/kernel_batched.dir/depend

