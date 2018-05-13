#
# Copyright 2016 WebAssembly Community Group participants
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

#Call with: cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=${EMSCRIPTEN}/cmake/Modules/Platform/Emscripten.cmake -G "Unix Makefiles" ../

cmake_minimum_required(VERSION 2.6)
project(ex1)

option(CODE_COVERAGE "Build with code coverage enabled" OFF)
option(WITH_EXCEPTIONS "Build with exceptions enabled" OFF)

if (${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
  set(COMPILER_IS_CLANG 1)
  set(COMPILER_IS_GNU 0)
  set(COMPILER_IS_MSVC 0)
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(COMPILER_IS_CLANG 0)
  set(COMPILER_IS_GNU 1)
  set(COMPILER_IS_MSVC 0)
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "MSVC")
  set(COMPILER_IS_CLANG 0)
  set(COMPILER_IS_GNU 0)
  set(COMPILER_IS_MSVC 1)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Emscripten")
  set(COMPILER_IS_CLANG 1)
  set(COMPILER_IS_GNU 0)
  set(COMPILER_IS_MSVC 0)
else ()
  set(COMPILER_IS_CLANG 0)
  set(COMPILER_IS_GNU 0)
  set(COMPILER_IS_MSVC 0)
endif ()

if (WIN32)
  check_symbol_exists(ENABLE_VIRTUAL_TERMINAL_PROCESSING "windows.h" HAVE_WIN32_VT100)
endif ()

if (EMSCRIPTEN)
  set(SIZEOF_SSIZE_T 4)
  set(SIZEOF_SIZE_T 4)
else ()
  include(CheckTypeSize)
  check_type_size(ssize_t SSIZE_T)
  check_type_size(size_t SIZEOF_SIZE_T)
endif ()

#configure_file(
#  ${WABT_SOURCE_DIR}/src/config.h.in
#  ${WABT_BINARY_DIR}/config.h
#)

include_directories(${WABT_SOURCE_DIR} ${WABT_BINARY_DIR})

if (COMPILER_IS_MSVC)
  # disable warning C4018: signed/unsigned mismatch
  # disable warning C4056, C4756: overflow in floating-point constant arithmetic
  #   seems to not like float compare w/ HUGE_VALF; bug?
  # disable warnings C4267 and C4244: conversion/truncation from larger to smaller type.
  # disable warning C4800: implicit conversion from larger int to bool
  add_definitions(-W3 -WX -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS)

  if (NOT WITH_EXCEPTIONS)
    # disable exception use in C++ library
    add_definitions(-D_HAS_EXCEPTIONS=0)
  endif ()
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  # disable -Wunused-parameter: this is really common when implementing
  #   interfaces, etc.
  # disable -Wpointer-arith: this is a GCC extension, and doesn't work in MSVC.
  add_definitions(
    -Wall -Wextra -Werror -Wno-unused-parameter -Wpointer-arith -g
    -Wold-style-cast -Wuninitialized
  )

  if (NOT WITH_EXCEPTIONS)
    add_definitions(-fno-exceptions)
  endif ()

  # Need to define __STDC_*_MACROS because C99 specifies that C++ shouldn't
  # define format (e.g. PRIu64) or limit (e.g. UINT32_MAX) macros without the
  # definition, and some libcs (e.g. glibc2.17 and earlier) follow that.
  add_definitions(-D__STDC_LIMIT_MACROS=1 -D__STDC_FORMAT_MACROS=1)

endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${WABT_SOURCE_DIR}/cmake)

add_custom_target(everything)

add_library(libglmatrix STATIC
  src/gl-matrix/dummy.c
  #src/gl-matrix/vec4.c
  #src/gl-matrix/mat4.c
  #src/gl-matrix/vec3.c
)
set_target_properties(libglmatrix PROPERTIES OUTPUT_NAME glmatrix)

if (NOT EMSCRIPTEN)
  function(wabt_executable name)
    # ARGV contains all arguments; remove the first one, ${name}, so it's just
    # a list of sources.
    list(REMOVE_AT ARGV 0)

    add_executable(${name} ${ARGV})
    add_dependencies(everything ${name})
    target_link_libraries(${name} libglmatrix)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD_REQUIRED ON)
    list(APPEND WABT_EXECUTABLES ${name})
    set(WABT_EXECUTABLES ${WABT_EXECUTABLES} PARENT_SCOPE)

#    add_custom_target(${name}-copy-to-bin ALL
#      COMMAND ${CMAKE_COMMAND} -E make_directory ${WABT_SOURCE_DIR}/bin
#      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${name}> ${WABT_SOURCE_DIR}/bin
#      DEPENDS ${name}
#    )
  endfunction()

  # ex1 executable
  wabt_executable(ex1 src/main.c)
  if (COMPILER_IS_CLANG OR COMPILER_IS_GNU)
    target_link_libraries(ex1 m)
  endif ()

  if (NOT CMAKE_VERSION VERSION_LESS "3.2")
    set(USES_TERMINAL USES_TERMINAL)
  endif ()

  # install
  install(TARGETS ${WABT_EXECUTABLES} DESTINATION bin)

else ()
  # emscripten stuff

  # just dump everything into one binary so we can reference it from JavaScript
  add_definitions(-Wno-warn-absolute-paths)
  add_executable(libex1js src/main.c)
  add_dependencies(everything libex1js)
  target_link_libraries(libex1js libglmatrix)
  set_target_properties(libex1js PROPERTIES OUTPUT_NAME libex1)

  set(LIBWABT_LINK_FLAGS
    --memory-init-file 0
    -s RESERVED_FUNCTION_POINTERS=10
    -s NO_EXIT_RUNTIME=1
    -s ALLOW_MEMORY_GROWTH=0
    -s ELIMINATE_DUPLICATE_FUNCTIONS=1
    -s WASM=1
    -Oz
    --llvm-lto 1
  )
  string(REPLACE ";" " " LIBWABT_LINK_FLAGS_STR "${LIBWABT_LINK_FLAGS}")

  set_target_properties(libex1js
    PROPERTIES
    LINK_FLAGS "${LIBWABT_LINK_FLAGS_STR}"
  )
endif ()