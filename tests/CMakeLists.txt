cmake_minimum_required(VERSION 3.16)
project(polyfit_tests LANGUAGES C CXX)

# -----------------------------------------------------------------------------
# Basic config
# -----------------------------------------------------------------------------
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable LTO if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)
if (LTO_SUPPORTED)
    message(STATUS "LTO supported; enabling interprocedural optimization")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
else ()
    message(WARNING "LTO not supported: ${LTO_ERROR}")
endif ()

# -----------------------------------------------------------------------------
# Fetch GoogleTest & Benchmark
# -----------------------------------------------------------------------------
CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION 1.17.0
        OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt OFF"
        "BUILD_GMOCK OFF"
)


CPMAddPackage(
        nanobench
        GITHUB_REPOSITORY martinus/nanobench
        GIT_TAG v4.3.11
)

# -----------------------------------------------------------------------------
# Option to turn on sanitizers
# -----------------------------------------------------------------------------
option(MONOFIT_ENABLE_SANITIZERS "Enable Address/UB sanitizers" ON)

# -----------------------------------------------------------------------------
# Helpers to add unit-test targets for multiple C++ standards
# -----------------------------------------------------------------------------
# Standards to test/build (C++17 and C++20). C++23 can be added later once mdspan handling is unified.
set(POLYFIT_TEST_STANDARDS 17 20)

function(add_polyfit_test_for_std BASENAME SRC STD)
    set(TGT "${BASENAME}_cxx${STD}")
    add_executable(${TGT} ${SRC})
    target_link_libraries(${TGT}
            PRIVATE polyfit::polyfit gtest_main
    )
    # Per-target standard selection
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD ${STD})
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${TGT} PROPERTY CXX_EXTENSIONS OFF)

    # common warnings & optimizations
    target_compile_options(${TGT} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic -O3 -march=native -ffast-math>
            $<$<CXX_COMPILER_ID:MSVC>:/W4 /O2 /arch:AVX2 /fp:fast>
    )

    if (MONOFIT_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TGT} PRIVATE -fsanitize=address,undefined)
        target_link_libraries(${TGT} PRIVATE -fsanitize=address,undefined)
    endif ()

    add_test(NAME ${TGT} COMMAND ${TGT})
endfunction()

function(add_polyfit_test_matrix NAME SRC)
    foreach(STD IN LISTS POLYFIT_TEST_STANDARDS)
        add_polyfit_test_for_std(${NAME} ${SRC} ${STD})
    endforeach()
endfunction()

# -----------------------------------------------------------------------------
# Helper to add a benchmark target
# -----------------------------------------------------------------------------
function(add_polyfit_benchmark NAME SRC)
    add_executable(${NAME} ${SRC})
    target_link_libraries(${NAME}
            PRIVATE
            polyfit::polyfit
            nanobench::nanobench
    )
    # Default to C++20 for the base target
    set_property(TARGET ${NAME} PROPERTY CXX_STANDARD 20)
    set_property(TARGET ${NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${NAME} PROPERTY CXX_EXTENSIONS OFF)

    target_compile_options(${NAME} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native -ffast-math>
            $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2>
    )
    # not a CTest test by default
endfunction()

# Benchmarks matrix across standards (where supported)
function(add_polyfit_benchmark_for_std BASENAME SRC STD)
    set(TGT "${BASENAME}_cxx${STD}")
    add_executable(${TGT} ${SRC})
    target_link_libraries(${TGT}
            PRIVATE
            polyfit::polyfit
            nanobench::nanobench
    )
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD ${STD})
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${TGT} PROPERTY CXX_EXTENSIONS OFF)
    target_compile_options(${TGT} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native -ffast-math>
            $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2>
    )
endfunction()

function(add_polyfit_benchmark_matrix NAME SRC)
    foreach(STD IN LISTS POLYFIT_TEST_STANDARDS)
        add_polyfit_benchmark_for_std(${NAME} ${SRC} ${STD})
    endforeach()
endfunction()

add_polyfit_test_matrix(test_horner test_horner.cpp)
add_polyfit_benchmark_matrix(bench_horner bench_horner.cpp)
add_polyfit_benchmark_matrix(bench_unroll_sweep bench_unroll_sweep.cpp)

add_polyfit_test_matrix(test1D test_1D.cpp)
add_polyfit_benchmark_matrix(bench1D bench_1D.cpp)


add_polyfit_test_matrix(testMany test_many.cpp)
add_polyfit_benchmark_matrix(benchMany bench_many.cpp)

add_polyfit_test_matrix(test_ND test_ND.cpp)
add_polyfit_benchmark_matrix(bench_ND bench_ND.cpp)

add_polyfit_test_matrix(test_general test_general.cpp)
add_polyfit_test_matrix(test_copy_move test_copy_move.cpp)
add_polyfit_test_matrix(test_scaling_alignment test_scaling_alignment.cpp)
add_polyfit_test_matrix(test_utils_bjorck_newton test_utils_bjorck_newton.cpp)
add_polyfit_test_matrix(test_utils_linspace test_utils_linspace.cpp)
add_polyfit_test_matrix(test_func_eval_many_padding test_func_eval_many_padding.cpp)
enable_testing()
