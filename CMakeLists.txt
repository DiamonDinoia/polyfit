cmake_minimum_required(VERSION 3.14)

project(polyfit
        VERSION 1.0.0
        DESCRIPTION "A fast polynomial fitting evaluation library"
        LANGUAGES CXX
)

# -----------------------------------------
# CPM (CMake Package Manager) bootstrap
# -----------------------------------------
set(CPM_DOWNLOAD_VERSION 0.42.0)
set(CPM_HASH_SUM "2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a")

if (CPM_SOURCE_CACHE)
    set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif (DEFINED ENV{CPM_SOURCE_CACHE})
    set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else ()
    set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif ()
get_filename_component(CPM_DOWNLOAD_LOCATION ${CPM_DOWNLOAD_LOCATION} ABSOLUTE)

file(DOWNLOAD
        https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
        ${CPM_DOWNLOAD_LOCATION}
        EXPECTED_HASH SHA256=${CPM_HASH_SUM}
        TLS_VERIFY ON
)
include(${CPM_DOWNLOAD_LOCATION})

# -----------------------------------------
# Dependencies (as subprojects so they can install with us)
# -----------------------------------------
CPMAddPackage(
        NAME xsimd
        GITHUB_REPOSITORY xtensor-stack/xsimd
        GIT_TAG 13.2.0
        OPTIONS
        "XSIMD_BUILD_TESTS OFF"
        "XSIMD_BUILD_BENCHMARKS OFF"
)

CPMAddPackage(
        NAME mdspan
        GITHUB_REPOSITORY kokkos/mdspan
        GIT_TAG mdspan-0.6.0
        OPTIONS
        "MDSPAN_ENABLE_TESTS OFF"
        "MDSPAN_ENABLE_EXAMPLES OFF"
        # Typically ON by default, but make it explicit:
        "MDSPAN_ENABLE_INSTALL ON"
)

CPMAddPackage(
        NAME poet
        GITHUB_REPOSITORY DiamonDinoia/poet
        GIT_TAG main
        OPTIONS
        "POET_BUILD_TESTS OFF"
        "BUILD_TESTING OFF"
)

# mdspan usually exports std::mdspan; if only mdspan::mdspan exists, create a compat alias.
if (NOT TARGET std::mdspan AND TARGET mdspan::mdspan)
    add_library(std::mdspan INTERFACE IMPORTED)
    target_link_libraries(std::mdspan INTERFACE mdspan::mdspan)
endif ()

# -----------------------------------------
# Library (header-only)
# -----------------------------------------
add_library(polyfit INTERFACE)
add_library(polyfit::polyfit ALIAS polyfit)

set_target_properties(polyfit PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

target_include_directories(polyfit INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# IMPORTANT:
# Only link deps in the BUILD interface; leave INSTALL interface empty.
# We'll re-attach deps for consumers inside polyfitConfig.cmake via find_dependency(...)
target_link_libraries(polyfit INTERFACE
        $<BUILD_INTERFACE:xsimd;std::mdspan;poet::poet>
)

# -----------------------------------------
# Install (headers + targets + package config)
# -----------------------------------------
include(GNUInstallDirs)

# Install public headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

# Install our interface target and export it
install(TARGETS polyfit
        EXPORT polyfitTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Export file that lists the installed targets
install(EXPORT polyfitTargets
        FILE polyfitTargets.cmake
        NAMESPACE polyfit::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyfit
)

# Configure the package config + version files
include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/polyfitConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/polyfitConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyfit
)

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/polyfitConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/polyfitConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/polyfitConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyfit
)

# -----------------------------------------
# Options: tests / examples
# -----------------------------------------
option(MONOFIT_BUILD_TESTS "Build tests" OFF)
if (MONOFIT_BUILD_TESTS)
    option(MONOFIT_ENABLE_SANITIZERS "Enable sanitizers for debugging (ASan/UBSan/etc.)" ON)
    enable_testing()
    add_subdirectory(tests)
endif ()

option(MONOFIT_BUILD_EXAMPLES "Build examples" ON)
if (MONOFIT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
