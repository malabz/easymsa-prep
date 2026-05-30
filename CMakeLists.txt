cmake_minimum_required(VERSION 3.20)
project(easymsa-prep VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type" FORCE)
endif()

find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package(LibLZMA REQUIRED)

set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_FUZZERS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(INSTALL_PROJECT OFF CACHE BOOL "" FORCE)
add_subdirectory(extern/miniz)

set(EASYMSA_PREP_MODULES
  src/cli.cpp
  src/config.cpp
  src/common.cpp
  src/sequence_utils.cpp
  src/input.cpp
  src/engine.cpp
  src/reports.cpp
)

add_library(easymsa_prep_core
  src/easymsa_prep.cpp
  src/xxhash.c
  ${EASYMSA_PREP_MODULES}
)

target_include_directories(easymsa_prep_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/extern/CLI11/include
    ${PROJECT_SOURCE_DIR}/extern/doctest/include
)

target_compile_definitions(easymsa_prep_core
  PUBLIC
    MINIZ_NO_ZLIB_COMPATIBLE_NAMES
    EASYMSA_PREP_VERSION="${PROJECT_VERSION}"
)

target_link_libraries(easymsa_prep_core
  PUBLIC
    miniz
    OpenMP::OpenMP_CXX
    ZLIB::ZLIB
    LibLZMA::LibLZMA
)

if(NOT MSVC)
  target_compile_options(easymsa_prep_core PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-pipe>
  )
endif()

add_executable(easymsa-prep src/main.cpp)
target_include_directories(easymsa-prep PRIVATE
  ${PROJECT_SOURCE_DIR}/extern/CLI11/include
)
target_link_libraries(easymsa-prep PRIVATE easymsa_prep_core)

install(TARGETS easymsa-prep RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
  add_executable(easymsa_prep_tests tests/test_core.cpp)
  target_include_directories(easymsa_prep_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/extern/doctest/include
  )
  target_link_libraries(easymsa_prep_tests PRIVATE easymsa_prep_core)
  add_test(NAME easymsa_prep_tests COMMAND easymsa_prep_tests)
endif()
