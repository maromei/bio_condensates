cmake_minimum_required(VERSION 3.13)
project(bio_condensates
  VERSION 1
  LANGUAGES CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

include(DuneMacros)

dune_project()
dune_enable_all_packages()

add_subdirectory(amdis)
add_subdirectory(src)

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
