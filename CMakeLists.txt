cmake_minimum_required(VERSION 3.14)

project(polyfit
        VERSION 1.0.0
        DESCRIPTION "A fast polynomial fitting evaluation library"
        LANGUAGES CXX)


# Download CPM (CMake Package Manager)
set(CPM_DOWNLOAD_VERSION 0.42.0)
set(CPM_HASH_SUM "2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a")

if (CPM_SOURCE_CACHE)
    set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif (DEFINED ENV{CPM_SOURCE_CACHE})
    set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else ()
    set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif ()

# Expand relative path. This is important if the provided path contains a tilde (~)
get_filename_component(CPM_DOWNLOAD_LOCATION ${CPM_DOWNLOAD_LOCATION} ABSOLUTE)

file(DOWNLOAD
        https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
        ${CPM_DOWNLOAD_LOCATION} EXPECTED_HASH SHA256=${CPM_HASH_SUM}
)

include(${CPM_DOWNLOAD_LOCATION})


# use CPM to fetch dependencies

CPMAddPackage(
        NAME xsimd
        GITHUB_REPOSITORY xtensor-stack/xsimd
        GIT_TAG 13.2.0
)

CPMAddPackage(
        mdspan
        GITHUB_REPOSITORY kokkos/mdspan
        GIT_TAG mdspan-0.6.0
)

# Create header-only library target
add_library(polyfit INTERFACE)
set_target_properties(polyfit PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
# Add alias for consistent naming
add_library(polyfit::polyfit ALIAS polyfit)
# Link dependencies to the polyfit library
target_link_libraries(polyfit INTERFACE xsimd std::mdspan)

# Set include directories for the library
target_include_directories(polyfit INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Install the library
include(GNUInstallDirs)

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

# Install the target
install(TARGETS polyfit
        EXPORT polyfitTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Create and install the config files
install(EXPORT polyfitTargets
        FILE polyfitTargets.cmake
        NAMESPACE polyfit::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyfit
)

# Generate config file
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

# Optional: Enable testing
option(MONOFIT_BUILD_TESTS "Build tests" ON)
if (MONOFIT_BUILD_TESTS)
    option(MONOFIT_ENABLE_SANITIZERS "Enable sanitizers for debugging (ASan/UBSan/etc.)" ON)
    enable_testing()
    add_subdirectory(tests)
endif ()

# Optional: Build examples
option(MONOFIT_BUILD_EXAMPLES "Build examples" ON)
if (MONOFIT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()