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
   ProgramOptions 
   sfml-graphics
   sfml-window
   sfml-system
   SimpleIni 
   tinyfiledialogs
   whereami)

add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND 
  ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_SOURCE_DIR}/resources ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (UNIX)
  set(INSTALL_BASEDIR "${CMAKE_INSTALL_PREFIX}/share/fmcomposer")

  install(PROGRAMS "$<TARGET_FILE:${EXECUTABLE_NAME}>" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" RENAME "fmcomposer.bin")
  install(FILES "${CMAKE_SOURCE_DIR}/fmcomposer.png"
          DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pixmaps")

  install(PROGRAMS "${CMAKE_SOURCE_DIR}/fmcomposer.sh" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" RENAME "fmcomposer")
  install(FILES "${CMAKE_SOURCE_DIR}/fmcomposer.desktop"
          DESTINATION "${CMAKE_INSTALL_PREFIX}/share/applications")

  # Update launch script/desktop file with path information
  install(CODE "file(APPEND \"${CMAKE_INSTALL_PREFIX}/bin/fmcomposer\" \"\nexec \\\"${CMAKE_INSTALL_PREFIX}/bin/fmcomposer.bin\\\" \\\"$@\\\" \\\"--appdir\\\" \\\"${INSTALL_BASEDIR}\\\"\n\")")
  install(CODE "file(APPEND \"${CMAKE_INSTALL_PREFIX}/share/applications/fmcomposer.desktop\" \"\nExec=\\\"${CMAKE_INSTALL_PREFIX}/bin/fmcomposer.bin\\\" \\\"$@\\\" \\\"--appdir\\\" \\\"${INSTALL_BASEDIR}\\\"\")")
  install(CODE "file(APPEND \"${CMAKE_INSTALL_PREFIX}/share/applications/fmcomposer.desktop\" \"\nIcon=\"${CMAKE_INSTALL_PREFIX}/share/pixmaps/fmcomposer.png\"\n\")")

  install(
    DIRECTORY "${CMAKE_SOURCE_DIR}/resources/instruments"
              "${CMAKE_SOURCE_DIR}/resources/songs"
              "${CMAKE_SOURCE_DIR}/resources/themes"
    DESTINATION ${INSTALL_BASEDIR}
    CONFIGURATIONS Release
  )
  install(
      FILES "${CMAKE_SOURCE_DIR}/resources/defaultkeyboards.ini"
            "${CMAKE_SOURCE_DIR}/resources/gmlist.ini"
            "${CMAKE_SOURCE_DIR}/resources/readme.txt"  
      DESTINATION "${INSTALL_BASEDIR}"
      CONFIGURATIONS Release
    )
endif()