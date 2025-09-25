cmake_minimum_required(VERSION 3.25)

project(FMComposer)
set(EXECUTABLE_NAME "FMComposer")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(MSVC)
  set(SFML2_DIR "${CMAKE_SOURCE_DIR}/libraries/SFML2")
  find_package(SFML2 REQUIRED)
else()
  find_package(SFML 2 COMPONENTS graphics window system REQUIRED)
  set(SFML2_LIBRARIES "sfml-graphics;sfml-window;sfml-system")
endif()

find_package(OpenGL REQUIRED)

if(MSVC)
    # Use static C runtime, means matching C runtime doesn't need to be on users box
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:Debug>")   
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /fp:fast")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:fast")
    # get the number of logical cores for parallel build
    cmake_host_system_information(RESULT LOGICAL_CORES QUERY NUMBER_OF_LOGICAL_CORES)
    math(EXPR COMPILE_CORES "${LOGICAL_CORES} - 1")  
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP${COMPILE_CORES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP${COMPILE_CORES}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS")
endif()

if (NOT MSVC AND ${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
endif ()

add_subdirectory(libraries)
include_directories(src)

if (MSVC)
  file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c" "src/fmcomposer.rc")
else()
  file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")
endif()

add_executable(${EXECUTABLE_NAME} ${SOURCES})

target_link_libraries(
   ${EXECUTABLE_NAME} 
   OpenGL::GL 
   portaudio_static
   portmidi
   ProgramOptions 
   ${SFML2_LIBRARIES}
   SimpleIni 
   tinyfiledialogs
   whereami)

add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin)

if (MSVC)
  target_include_directories(${EXECUTABLE_NAME} SYSTEM PRIVATE ${SFML2_INCLUDE_DIRS})
  add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND_EXPAND_LISTS COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SFML2_DLLS} ${CMAKE_BINARY_DIR}/bin)
endif()

add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:${EXECUTABLE_NAME}>" ${CMAKE_BINARY_DIR}/bin)
add_custom_command( TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND 
  ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/bin)

if (UNIX)
  set(INSTALL_BASEDIR "${CMAKE_INSTALL_PREFIX}/share/fmcomposer")

  install(PROGRAMS "$<TARGET_FILE:${EXECUTABLE_NAME}>" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" RENAME "fmcomposer.bin")
  install(FILES "${CMAKE_SOURCE_DIR}/src/fmcomposer.png"
          DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pixmaps")

  install(PROGRAMS "${CMAKE_SOURCE_DIR}/src/fmcomposer.sh" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" RENAME "fmcomposer")
  install(FILES "${CMAKE_SOURCE_DIR}/src/fmcomposer.desktop"
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