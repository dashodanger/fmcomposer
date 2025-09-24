cmake_minimum_required(VERSION 3.25)

if(NOT DEFINED CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build (Debug or Release)")
endif()

project(FMComposer)
set(EXECUTABLE_NAME "FMComposer")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SFML 2 COMPONENTS graphics window system REQUIRED)
find_package(OpenGL REQUIRED)

if (${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
endif ()

add_subdirectory(libraries)
include_directories(src)

file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(${EXECUTABLE_NAME} ${SOURCES})

target_link_libraries(
   ${EXECUTABLE_NAME} 
   OpenGL::GL 
   portaudio_static
   portmidi 
   sfml-graphics
   sfml-window
   sfml-system
   SimpleIni 
   tinyfiledialogs)

add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND 
  ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_SOURCE_DIR}/resources ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})