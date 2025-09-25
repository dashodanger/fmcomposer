set(SFML2_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")

if (${CMAKE_SIZEOF_VOID_P} MATCHES 8)
  set(SFML2_ARCH "x64")
else()
  set(SFML2_ARCH "x86")
endif()

set(SFML2_LIBRARIES "${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/lib/freetype.lib;${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/lib/sfml-graphics.lib;${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/lib/sfml-system.lib;${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/lib/sfml-window.lib;${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/lib/sfml-main.lib")
set(SFML2_DLLS "")
list (APPEND SFML2_DLLS "${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/bin/sfml-graphics-2.dll" "${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/bin/sfml-system-2.dll" "${CMAKE_CURRENT_LIST_DIR}/${SFML2_ARCH}/bin/sfml-window-2.dll")

string(STRIP "${SFML2_LIBRARIES}" SFML2_LIBRARIES)