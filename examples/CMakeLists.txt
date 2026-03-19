cmake_minimum_required(VERSION 3.14)
project(polyfit_examples VERSION 1.0.0 LANGUAGES CXX)

include(PolyfitWarnings)
include(PolyfitSanitizers)
include(PolyfitStaticAnalysis)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(POLYFIT_EXAMPLE_STANDARDS 17)
foreach(_std IN ITEMS 20 23)
    check_cxx_compiler_flag("-std=c++${_std}" POLYFIT_EXAMPLE_HAS_CXX${_std})
    if(POLYFIT_EXAMPLE_HAS_CXX${_std})
        list(APPEND POLYFIT_EXAMPLE_STANDARDS ${_std})
    endif()
endforeach()

if("20" IN_LIST POLYFIT_EXAMPLE_STANDARDS)
    set(COMPILER_SUPPORTS_CXX20 TRUE)
else()
    set(COMPILER_SUPPORTS_CXX20 FALSE)
    message(STATUS "Compiler does not support C++20; compile-time examples will be skipped")
endif()

set(SOURCES
        runtime_1D_fixed.cpp
        runtime_1D_eps.cpp
        runtime_ND_fixed.cpp
        pack.cpp
        tags.cpp
)
set(SOURCES_CXX20)
if(COMPILER_SUPPORTS_CXX20)
    list(APPEND SOURCES_CXX20
            compiletime_1D_fixed.cpp
            compiletime_1D_eps.cpp
    )
endif()

function(polyfit_configure_example_target target)
    target_link_libraries(${target} PRIVATE polyfit::polyfit)
    if (POLYFIT_STRICT_WARNINGS)
        polyfit_enable_warnings(${target})
    endif ()
    polyfit_enable_sanitizers(${target})
    polyfit_configure_static_analysis(${target})
    if (POLYFIT_ENABLE_NATIVE_TUNING)
        target_compile_options(${target} PRIVATE
                $<$<CXX_COMPILER_ID:GNU,Clang>:-march=native>
        )
    endif ()
endfunction()

set(EXAMPLE_TARGETS "")
foreach(src IN LISTS SOURCES SOURCES_CXX20)
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    polyfit_configure_example_target(${name})
    if(src IN_LIST SOURCES_CXX20)
        set_property(TARGET ${name} PROPERTY CXX_STANDARD 20)
        set_property(TARGET ${name} PROPERTY CXX_STANDARD_REQUIRED ON)
    endif()
    list(APPEND EXAMPLE_TARGETS ${name})
endforeach()

foreach(std IN LISTS POLYFIT_EXAMPLE_STANDARDS)
    set(name portable_trig_cxx${std})
    add_executable(${name} portable_trig.cpp)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD ${std})
    set_property(TARGET ${name} PROPERTY CXX_STANDARD_REQUIRED ON)
    polyfit_configure_example_target(${name})
    list(APPEND EXAMPLE_TARGETS ${name})
endforeach()

set(RUN_COMMANDS
        COMMAND ${CMAKE_COMMAND} -E echo "=== Building and running all examples ==="
)
foreach(exe IN LISTS EXAMPLE_TARGETS)
    list(APPEND RUN_COMMANDS
            COMMAND ${CMAKE_COMMAND} -E echo "--> Running ${exe}"
            COMMAND $<TARGET_FILE:${exe}>
    )
endforeach()

add_custom_target(run_examples
        DEPENDS ${EXAMPLE_TARGETS}
        USES_TERMINAL
        COMMENT "Building and running all examples"
        ${RUN_COMMANDS}
)
