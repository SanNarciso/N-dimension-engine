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
include math_lib/CMakeFiles/math_lib.dir/depend.make

# Include the progress variables for this target.
include math_lib/CMakeFiles/math_lib.dir/progress.make

# Include the compile flags for this target's objects.
include math_lib/CMakeFiles/math_lib.dir/flags.make

math_lib/CMakeFiles/math_lib.dir/__/base.obj: math_lib/CMakeFiles/math_lib.dir/flags.make
math_lib/CMakeFiles/math_lib.dir/__/base.obj: math_lib/CMakeFiles/math_lib.dir/includes_CXX.rsp
math_lib/CMakeFiles/math_lib.dir/__/base.obj: ../base.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object math_lib/CMakeFiles/math_lib.dir/__/base.obj"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\math_lib.dir\__\base.obj -c "C:\University\IV\Computer graphics project\base.cpp"

math_lib/CMakeFiles/math_lib.dir/__/base.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/math_lib.dir/__/base.i"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\University\IV\Computer graphics project\base.cpp" > CMakeFiles\math_lib.dir\__\base.i

math_lib/CMakeFiles/math_lib.dir/__/base.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/math_lib.dir/__/base.s"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\University\IV\Computer graphics project\base.cpp" -o CMakeFiles\math_lib.dir\__\base.s

# Object files for target math_lib
math_lib_OBJECTS = \
"CMakeFiles/math_lib.dir/__/base.obj"

# External object files for target math_lib
math_lib_EXTERNAL_OBJECTS =

math_lib/libmath_lib.a: math_lib/CMakeFiles/math_lib.dir/__/base.obj
math_lib/libmath_lib.a: math_lib/CMakeFiles/math_lib.dir/build.make
math_lib/libmath_lib.a: math_lib/CMakeFiles/math_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libmath_lib.a"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && $(CMAKE_COMMAND) -P CMakeFiles\math_lib.dir\cmake_clean_target.cmake
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\math_lib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
math_lib/CMakeFiles/math_lib.dir/build: math_lib/libmath_lib.a

.PHONY : math_lib/CMakeFiles/math_lib.dir/build

math_lib/CMakeFiles/math_lib.dir/clean:
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\math_lib && $(CMAKE_COMMAND) -P CMakeFiles\math_lib.dir\cmake_clean.cmake
.PHONY : math_lib/CMakeFiles/math_lib.dir/clean

math_lib/CMakeFiles/math_lib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\University\IV\Computer graphics project" "C:\University\IV\Computer graphics project\math_lib" "C:\University\IV\Computer graphics project\cmake-build-debug" "C:\University\IV\Computer graphics project\cmake-build-debug\math_lib" "C:\University\IV\Computer graphics project\cmake-build-debug\math_lib\CMakeFiles\math_lib.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : math_lib/CMakeFiles/math_lib.dir/depend

