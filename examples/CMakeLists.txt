cmake_minimum_required(VERSION 3.14)
project(polyfit_examples VERSION 1.0.0 LANGUAGES CXX)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++20" COMPILER_SUPPORTS_CXX20)
if (COMPILER_SUPPORTS_CXX20)
    message(STATUS "Enabling C++20")
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(WARNING "Compiler does NOT support C++20; using default C++ standard")
endif()

# 1. List your source files
set(SOURCES
        runtime_1D_fixed.cpp
        runtime_1D_eps.cpp
        compiletime_1D_fixed.cpp
        runtime_ND_fixed.cpp
        runtime_ND_eps.cpp
)
if (COMPILER_SUPPORTS_CXX20)
    list(APPEND SOURCES compiletime_1D_eps.cpp)
endif()

# 2. Generate one executable target per source, collect their names:
set(EXAMPLE_TARGETS "")
foreach(src IN LISTS SOURCES)
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE polyfit::polyfit)
    list(APPEND EXAMPLE_TARGETS ${name})
endforeach()

# 3. Build up all of the “run” commands in a list:
set(RUN_COMMANDS
        # This first COMMAND makes USES_TERMINAL happy:
        COMMAND ${CMAKE_COMMAND} -E echo "=== Building and running all examples ==="
)
foreach(exe IN LISTS EXAMPLE_TARGETS)
    list(APPEND RUN_COMMANDS
            COMMAND ${CMAKE_COMMAND} -E echo "--> Running ${exe}"
            COMMAND $<TARGET_FILE:${exe}>
    )
endforeach()

# 4. Single custom target that (a) depends on all the examples, then (b) runs them
add_custom_target(run_examples
        DEPENDS ${EXAMPLE_TARGETS}
        USES_TERMINAL
        COMMENT "Building and running all examples"
        ${RUN_COMMANDS}
)
