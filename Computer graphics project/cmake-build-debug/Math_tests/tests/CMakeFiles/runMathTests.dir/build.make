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
include Math_tests/tests/CMakeFiles/runMathTests.dir/depend.make

# Include the progress variables for this target.
include Math_tests/tests/CMakeFiles/runMathTests.dir/progress.make

# Include the compile flags for this target's objects.
include Math_tests/tests/CMakeFiles/runMathTests.dir/flags.make

Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.obj: Math_tests/tests/CMakeFiles/runMathTests.dir/flags.make
Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.obj: Math_tests/tests/CMakeFiles/runMathTests.dir/includes_CXX.rsp
Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.obj: ../Math_tests/tests/testing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.obj"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\tests && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\runMathTests.dir\testing.cpp.obj -c "C:\University\IV\Computer graphics project\Math_tests\tests\testing.cpp"

Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMathTests.dir/testing.cpp.i"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\tests && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\University\IV\Computer graphics project\Math_tests\tests\testing.cpp" > CMakeFiles\runMathTests.dir\testing.cpp.i

Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMathTests.dir/testing.cpp.s"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\tests && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\University\IV\Computer graphics project\Math_tests\tests\testing.cpp" -o CMakeFiles\runMathTests.dir\testing.cpp.s

# Object files for target runMathTests
runMathTests_OBJECTS = \
"CMakeFiles/runMathTests.dir/testing.cpp.obj"

# External object files for target runMathTests
runMathTests_EXTERNAL_OBJECTS =

Math_tests/tests/runMathTests.exe: Math_tests/tests/CMakeFiles/runMathTests.dir/testing.cpp.obj
Math_tests/tests/runMathTests.exe: Math_tests/tests/CMakeFiles/runMathTests.dir/build.make
Math_tests/tests/runMathTests.exe: lib/libgtest.a
Math_tests/tests/runMathTests.exe: lib/libgtest_main.a
Math_tests/tests/runMathTests.exe: lib/libgtest.a
Math_tests/tests/runMathTests.exe: Math_tests/tests/CMakeFiles/runMathTests.dir/linklibs.rsp
Math_tests/tests/runMathTests.exe: Math_tests/tests/CMakeFiles/runMathTests.dir/objects1.rsp
Math_tests/tests/runMathTests.exe: Math_tests/tests/CMakeFiles/runMathTests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable runMathTests.exe"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\runMathTests.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Math_tests/tests/CMakeFiles/runMathTests.dir/build: Math_tests/tests/runMathTests.exe

.PHONY : Math_tests/tests/CMakeFiles/runMathTests.dir/build

Math_tests/tests/CMakeFiles/runMathTests.dir/clean:
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\tests && $(CMAKE_COMMAND) -P CMakeFiles\runMathTests.dir\cmake_clean.cmake
.PHONY : Math_tests/tests/CMakeFiles/runMathTests.dir/clean

Math_tests/tests/CMakeFiles/runMathTests.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\University\IV\Computer graphics project" "C:\University\IV\Computer graphics project\Math_tests\tests" "C:\University\IV\Computer graphics project\cmake-build-debug" "C:\University\IV\Computer graphics project\cmake-build-debug\Math_tests\tests" "C:\University\IV\Computer graphics project\cmake-build-debug\Math_tests\tests\CMakeFiles\runMathTests.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : Math_tests/tests/CMakeFiles/runMathTests.dir/depend

