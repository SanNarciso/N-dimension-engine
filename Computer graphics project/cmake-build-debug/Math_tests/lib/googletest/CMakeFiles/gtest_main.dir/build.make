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
include Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/depend.make

# Include the progress variables for this target.
include Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/progress.make

# Include the compile flags for this target's objects.
include Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/flags.make

Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj: Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/flags.make
Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj: Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/includes_CXX.rsp
Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj: ../Math_tests/lib/googletest/src/gtest_main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\gtest_main.dir\src\gtest_main.cc.obj -c "C:\University\IV\Computer graphics project\Math_tests\lib\googletest\src\gtest_main.cc"

Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gtest_main.dir/src/gtest_main.cc.i"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\University\IV\Computer graphics project\Math_tests\lib\googletest\src\gtest_main.cc" > CMakeFiles\gtest_main.dir\src\gtest_main.cc.i

Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gtest_main.dir/src/gtest_main.cc.s"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && D:\MinGW\mingw32\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\University\IV\Computer graphics project\Math_tests\lib\googletest\src\gtest_main.cc" -o CMakeFiles\gtest_main.dir\src\gtest_main.cc.s

# Object files for target gtest_main
gtest_main_OBJECTS = \
"CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj"

# External object files for target gtest_main
gtest_main_EXTERNAL_OBJECTS =

lib/libgtest_main.a: Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/src/gtest_main.cc.obj
lib/libgtest_main.a: Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/build.make
lib/libgtest_main.a: Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\University\IV\Computer graphics project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\..\..\lib\libgtest_main.a"
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && $(CMAKE_COMMAND) -P CMakeFiles\gtest_main.dir\cmake_clean_target.cmake
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\gtest_main.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/build: lib/libgtest_main.a

.PHONY : Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/build

Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/clean:
	cd /d C:\UNIVER~1\IV\COMPUT~1\CMAKE-~1\MATH_T~1\lib\GOOGLE~2 && $(CMAKE_COMMAND) -P CMakeFiles\gtest_main.dir\cmake_clean.cmake
.PHONY : Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/clean

Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\University\IV\Computer graphics project" "C:\University\IV\Computer graphics project\Math_tests\lib\googletest" "C:\University\IV\Computer graphics project\cmake-build-debug" "C:\University\IV\Computer graphics project\cmake-build-debug\Math_tests\lib\googletest" "C:\University\IV\Computer graphics project\cmake-build-debug\Math_tests\lib\googletest\CMakeFiles\gtest_main.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : Math_tests/lib/googletest/CMakeFiles/gtest_main.dir/depend

