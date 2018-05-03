# Install script for directory: /home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ASAGI

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/asagi_nompi.pc")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/usr/local/lib:/lrz/sys/libraries/netcdf/4.3.3/intel/serial_160/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/libasagi_nompi.so.1.0"
    "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/libasagi_nompi.so.1"
    "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/libasagi_nompi.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libasagi_nompi.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/lrz/sys/libraries/netcdf/4.3.3/intel/serial_160/lib:::::::::::::::"
           NEW_RPATH "/usr/local/lib:/lrz/sys/libraries/netcdf/4.3.3/intel/serial_160/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ASAGI/include/asagi.h"
    "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/ASAGI/include/asagi.f90"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/allocator/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/cache/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/fortran/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/grid/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/io/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/perf/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/magic/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/numa/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/threads/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/transfer/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/types/cmake_install.cmake")
  include("/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/documentation/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/hpc/pr63so/di29zux/ExaHyPE-Engine/3d/easi/libs/asagi_build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
