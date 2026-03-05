cmake_minimum_required(VERSION 3.14)

project(polyfit
        VERSION 1.0.0
        DESCRIPTION "A fast polynomial fitting evaluation library"
        LANGUAGES CXX
)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CTest)

# Backward-compatible aliases for legacy option names.
if (DEFINED MONOFIT_BUILD_TESTS AND NOT DEFINED POLYFIT_BUILD_TESTS)
    message(DEPRECATION "MONOFIT_BUILD_TESTS is deprecated; use POLYFIT_BUILD_TESTS")
    set(POLYFIT_BUILD_TESTS "${MONOFIT_BUILD_TESTS}" CACHE BOOL "Build tests" FORCE)
endif ()
if (DEFINED MONOFIT_BUILD_EXAMPLES AND NOT DEFINED POLYFIT_BUILD_EXAMPLES)
    message(DEPRECATION "MONOFIT_BUILD_EXAMPLES is deprecated; use POLYFIT_BUILD_EXAMPLES")
    set(POLYFIT_BUILD_EXAMPLES "${MONOFIT_BUILD_EXAMPLES}" CACHE BOOL "Build examples" FORCE)
endif ()
if (DEFINED MONOFIT_ENABLE_SANITIZERS AND NOT DEFINED POLYFIT_ENABLE_SANITIZERS)
    message(DEPRECATION "MONOFIT_ENABLE_SANITIZERS is deprecated; use POLYFIT_ENABLE_SANITIZERS")
    set(POLYFIT_ENABLE_SANITIZERS "${MONOFIT_ENABLE_SANITIZERS}" CACHE BOOL "Enable sanitizers" FORCE)
endif ()

if (NOT DEFINED POLYFIT_BUILD_TESTS)
    set(POLYFIT_BUILD_TESTS "${BUILD_TESTING}")
endif ()

option(POLYFIT_BUILD_TESTS "Build tests" "${POLYFIT_BUILD_TESTS}")
option(POLYFIT_BUILD_EXAMPLES "Build examples" ON)
option(POLYFIT_STRICT_WARNINGS "Enable strict warning profile for project targets" ON)
option(POLYFIT_ENABLE_NATIVE_TUNING "Enable -march=native tuning for benchmark/example targets" OFF)
option(POLYFIT_ENABLE_COVERAGE "Enable coverage instrumentation for project/test targets" OFF)

include(PolyfitWarnings)
include(PolyfitSanitizers)
include(PolyfitStaticAnalysis)
include(PolyfitDocs)
include(PolyfitAmalgamate)

# -----------------------------------------
# CPM (CMake Package Manager) bootstrap
# -----------------------------------------
set(CPM_DOWNLOAD_VERSION 0.42.1)
set(CPM_HASH_SUM "f3a6dcc6a04ce9e7f51a127307fa4f699fb2bade357a8eb4c5b45df76e1dc6a5")

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
        GIT_TAG 14.0.0
        OPTIONS
        "XSIMD_BUILD_TESTS OFF"
        "XSIMD_BUILD_BENCHMARKS OFF"
)

# Only download kokkos/mdspan when the compiler lacks native std::mdspan.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #include <version>
    #if !defined(__cpp_lib_mdspan) || __cpp_lib_mdspan < 202207L
    #error no native mdspan
    #endif
    int main() {}
" POLYFIT_HAS_NATIVE_MDSPAN)

if (NOT POLYFIT_HAS_NATIVE_MDSPAN)
    CPMAddPackage(
            NAME mdspan
            GITHUB_REPOSITORY kokkos/mdspan
            GIT_TAG mdspan-0.6.0
            OPTIONS
            "MDSPAN_ENABLE_TESTS OFF"
            "MDSPAN_ENABLE_EXAMPLES OFF"
            "MDSPAN_ENABLE_INSTALL ON"
    )
endif ()

CPMAddPackage(
        NAME poet
        GITHUB_REPOSITORY DiamonDinoia/poet
        GIT_TAG dfb54db
        OPTIONS
        "POET_BUILD_TESTS OFF"
        "BUILD_TESTING OFF"
)

# When using the kokkos CPM package: create std::mdspan alias if needed.
# When using native mdspan: create a no-op std::mdspan target so link lines stay uniform.
if (NOT TARGET std::mdspan)
    if (TARGET mdspan::mdspan)
        add_library(std::mdspan INTERFACE IMPORTED)
        target_link_libraries(std::mdspan INTERFACE mdspan::mdspan)
    else ()
        add_library(std::mdspan INTERFACE IMPORTED)
    endif ()
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

# Keep third-party headers as SYSTEM for consumers so strict warning profiles
# apply to project code rather than dependency internals.
foreach(_dep_target IN ITEMS xsimd std::mdspan poet::poet)
    if (TARGET ${_dep_target})
        get_target_property(_dep_includes ${_dep_target} INTERFACE_INCLUDE_DIRECTORIES)
        if (_dep_includes)
            target_include_directories(polyfit SYSTEM INTERFACE "${_dep_includes}")
        endif ()
    endif ()
endforeach ()

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
if (POLYFIT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

if (POLYFIT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
