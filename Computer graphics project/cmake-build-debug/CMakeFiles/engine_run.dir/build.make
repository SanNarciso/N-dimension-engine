# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.19

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\CLion 2021.1.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\CLion 2021.1.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\University\IV\Computer graphics project"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\University\IV\Computer graphics project\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/Engine_run.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Engine_run.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Engine_run.dir/flags.make

CMakeFiles/Engine_run.dir/main.cpp.obj: CMakeFiles/Engine_run.dir/flags.make
CMakeFiles/Engine_run.dir/main.cpp.obj: CMakeFiles/Engine_run.dir/includes_CXX.rsp
CMakeFiles/Engine_run.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Engine_run.dir/main.cpp.obj"
	D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Engine_run.dir\main.cpp.obj -c "C:\University\IV\Computer graphics project\main.cpp"

CMakeFiles/Engine_run.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Engine_run.dir/main.cpp.i"
	D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\University\IV\Computer graphics project\main.cpp" > CMakeFiles\Engine_run.dir\main.cpp.i

CMakeFiles/Engine_run.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Engine_run.dir/main.cpp.s"
	D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\University\IV\Computer graphics project\main.cpp" -o CMakeFiles\Engine_run.dir\main.cpp.s

# Object files for target Engine_run
Engine_run_OBJECTS = \
"CMakeFiles/Engine_run.dir/main.cpp.obj"

# External object files for target Engine_run
Engine_run_EXTERNAL_OBJECTS =

Engine_run.exe: CMakeFiles/Engine_run.dir/main.cpp.obj
Engine_run.exe: CMakeFiles/Engine_run.dir/build.make
Engine_run.exe: Engine/libEngine_lib.a
Engine_run.exe: CMakeFiles/Engine_run.dir/linklibs.rsp
Engine_run.exe: CMakeFiles/Engine_run.dir/objects1.rsp
Engine_run.exe: CMakeFiles/Engine_run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Engine_run.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Engine_run.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Engine_run.dir/build: Engine_run.exe

.PHONY : CMakeFiles/Engine_run.dir/build

CMakeFiles/Engine_run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Engine_run.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Engine_run.dir/clean

CMakeFiles/Engine_run.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\University\IV\Computer graphics project" "C:\University\IV\Computer graphics project" "C:\University\IV\Computer graphics project\cmake-build-debug" "C:\University\IV\Computer graphics project\cmake-build-debug" "C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles\engine_run.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/Engine_run.dir/depend

