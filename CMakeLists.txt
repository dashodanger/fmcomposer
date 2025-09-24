cmake_minimum_required(VERSION 3.25)

# Set default build type to Release
if(NOT DEFINED CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build (Debug or Release)")
endif()

project(FMComposer)

# Detect and enable C++11 support
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SFML 2 COMPONENTS graphics window system REQUIRED)
find_package(OpenGL REQUIRED)

add_subdirectory(libraries)
include_directories(src)

# Add source files
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")

# Set executable output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Set executable name
set(EXECUTABLE_NAME "FMComposer")
add_executable(${EXECUTABLE_NAME} ${SOURCES})

target_link_libraries(${EXECUTABLE_NAME} OpenGL::GL portaudio portmidi sfml-graphics sfml-window sfml-system SimpleIni tinyfiledialogs)