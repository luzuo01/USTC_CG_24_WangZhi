# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wangzhi/Desktop/USTC_CG_24/Framework2D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wangzhi/Desktop/USTC_CG_24/Framework2D

# Include any dependencies generated for this target.
include src/view/CMakeFiles/view.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/view/CMakeFiles/view.dir/compiler_depend.make

# Include the progress variables for this target.
include src/view/CMakeFiles/view.dir/progress.make

# Include the compile flags for this target's objects.
include src/view/CMakeFiles/view.dir/flags.make

src/view/CMakeFiles/view.dir/comp_canvas.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/comp_canvas.cpp.o: src/view/comp_canvas.cpp
src/view/CMakeFiles/view.dir/comp_canvas.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/view/CMakeFiles/view.dir/comp_canvas.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/comp_canvas.cpp.o -MF CMakeFiles/view.dir/comp_canvas.cpp.o.d -o CMakeFiles/view.dir/comp_canvas.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_canvas.cpp

src/view/CMakeFiles/view.dir/comp_canvas.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/comp_canvas.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_canvas.cpp > CMakeFiles/view.dir/comp_canvas.cpp.i

src/view/CMakeFiles/view.dir/comp_canvas.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/comp_canvas.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_canvas.cpp -o CMakeFiles/view.dir/comp_canvas.cpp.s

src/view/CMakeFiles/view.dir/comp_image.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/comp_image.cpp.o: src/view/comp_image.cpp
src/view/CMakeFiles/view.dir/comp_image.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/view/CMakeFiles/view.dir/comp_image.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/comp_image.cpp.o -MF CMakeFiles/view.dir/comp_image.cpp.o.d -o CMakeFiles/view.dir/comp_image.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_image.cpp

src/view/CMakeFiles/view.dir/comp_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/comp_image.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_image.cpp > CMakeFiles/view.dir/comp_image.cpp.i

src/view/CMakeFiles/view.dir/comp_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/comp_image.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/comp_image.cpp -o CMakeFiles/view.dir/comp_image.cpp.s

src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o: src/view/shapes/ellipse.cpp
src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o -MF CMakeFiles/view.dir/shapes/ellipse.cpp.o.d -o CMakeFiles/view.dir/shapes/ellipse.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/ellipse.cpp

src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/shapes/ellipse.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/ellipse.cpp > CMakeFiles/view.dir/shapes/ellipse.cpp.i

src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/shapes/ellipse.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/ellipse.cpp -o CMakeFiles/view.dir/shapes/ellipse.cpp.s

src/view/CMakeFiles/view.dir/shapes/line.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/shapes/line.cpp.o: src/view/shapes/line.cpp
src/view/CMakeFiles/view.dir/shapes/line.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/view/CMakeFiles/view.dir/shapes/line.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/shapes/line.cpp.o -MF CMakeFiles/view.dir/shapes/line.cpp.o.d -o CMakeFiles/view.dir/shapes/line.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/line.cpp

src/view/CMakeFiles/view.dir/shapes/line.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/shapes/line.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/line.cpp > CMakeFiles/view.dir/shapes/line.cpp.i

src/view/CMakeFiles/view.dir/shapes/line.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/shapes/line.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/line.cpp -o CMakeFiles/view.dir/shapes/line.cpp.s

src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o: src/view/shapes/polygon.cpp
src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o -MF CMakeFiles/view.dir/shapes/polygon.cpp.o.d -o CMakeFiles/view.dir/shapes/polygon.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/polygon.cpp

src/view/CMakeFiles/view.dir/shapes/polygon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/shapes/polygon.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/polygon.cpp > CMakeFiles/view.dir/shapes/polygon.cpp.i

src/view/CMakeFiles/view.dir/shapes/polygon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/shapes/polygon.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/polygon.cpp -o CMakeFiles/view.dir/shapes/polygon.cpp.s

src/view/CMakeFiles/view.dir/shapes/rect.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/shapes/rect.cpp.o: src/view/shapes/rect.cpp
src/view/CMakeFiles/view.dir/shapes/rect.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/view/CMakeFiles/view.dir/shapes/rect.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/shapes/rect.cpp.o -MF CMakeFiles/view.dir/shapes/rect.cpp.o.d -o CMakeFiles/view.dir/shapes/rect.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/rect.cpp

src/view/CMakeFiles/view.dir/shapes/rect.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/shapes/rect.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/rect.cpp > CMakeFiles/view.dir/shapes/rect.cpp.i

src/view/CMakeFiles/view.dir/shapes/rect.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/shapes/rect.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/shapes/rect.cpp -o CMakeFiles/view.dir/shapes/rect.cpp.s

src/view/CMakeFiles/view.dir/window.cpp.o: src/view/CMakeFiles/view.dir/flags.make
src/view/CMakeFiles/view.dir/window.cpp.o: src/view/window.cpp
src/view/CMakeFiles/view.dir/window.cpp.o: src/view/CMakeFiles/view.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/view/CMakeFiles/view.dir/window.cpp.o"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/view/CMakeFiles/view.dir/window.cpp.o -MF CMakeFiles/view.dir/window.cpp.o.d -o CMakeFiles/view.dir/window.cpp.o -c /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/window.cpp

src/view/CMakeFiles/view.dir/window.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/view.dir/window.cpp.i"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/window.cpp > CMakeFiles/view.dir/window.cpp.i

src/view/CMakeFiles/view.dir/window.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/view.dir/window.cpp.s"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/window.cpp -o CMakeFiles/view.dir/window.cpp.s

# Object files for target view
view_OBJECTS = \
"CMakeFiles/view.dir/comp_canvas.cpp.o" \
"CMakeFiles/view.dir/comp_image.cpp.o" \
"CMakeFiles/view.dir/shapes/ellipse.cpp.o" \
"CMakeFiles/view.dir/shapes/line.cpp.o" \
"CMakeFiles/view.dir/shapes/polygon.cpp.o" \
"CMakeFiles/view.dir/shapes/rect.cpp.o" \
"CMakeFiles/view.dir/window.cpp.o"

# External object files for target view
view_EXTERNAL_OBJECTS =

libs/libview_d.a: src/view/CMakeFiles/view.dir/comp_canvas.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/comp_image.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/shapes/ellipse.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/shapes/line.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/shapes/polygon.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/shapes/rect.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/window.cpp.o
libs/libview_d.a: src/view/CMakeFiles/view.dir/build.make
libs/libview_d.a: src/view/CMakeFiles/view.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/wangzhi/Desktop/USTC_CG_24/Framework2D/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX static library ../../libs/libview_d.a"
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && $(CMAKE_COMMAND) -P CMakeFiles/view.dir/cmake_clean_target.cmake
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/view.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/view/CMakeFiles/view.dir/build: libs/libview_d.a
.PHONY : src/view/CMakeFiles/view.dir/build

src/view/CMakeFiles/view.dir/clean:
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view && $(CMAKE_COMMAND) -P CMakeFiles/view.dir/cmake_clean.cmake
.PHONY : src/view/CMakeFiles/view.dir/clean

src/view/CMakeFiles/view.dir/depend:
	cd /Users/wangzhi/Desktop/USTC_CG_24/Framework2D && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wangzhi/Desktop/USTC_CG_24/Framework2D /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view /Users/wangzhi/Desktop/USTC_CG_24/Framework2D /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view /Users/wangzhi/Desktop/USTC_CG_24/Framework2D/src/view/CMakeFiles/view.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/view/CMakeFiles/view.dir/depend

