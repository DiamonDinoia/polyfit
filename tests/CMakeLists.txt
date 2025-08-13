cmake_minimum_required(VERSION 3.16)
project(polyfit_tests LANGUAGES C CXX)

# -----------------------------------------------------------------------------
# Basic config
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
        NAME benchmark
        GITHUB_REPOSITORY google/benchmark
        VERSION 1.9.4
        OPTIONS
        "BENCHMARK_ENABLE_TESTING OFF"
        "BENCHMARK_ENABLE_INSTALL OFF"
        "BENCHMARK_ENABLE_GTEST_TESTS OFF"
        "BENCHMARK_ENABLE_LTO ON"
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
# Helper to add a unit‐test target
# -----------------------------------------------------------------------------
function(add_polyfit_test NAME SRC)
    add_executable(${NAME} ${SRC})
    target_link_libraries(${NAME}
            PRIVATE polyfit::polyfit gtest_main
    )
    target_compile_features(${NAME} PRIVATE cxx_std_20)


    # common warnings & optimizations
    target_compile_options(${NAME} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic -O3 -march=native -ffast-math>
            $<$<CXX_COMPILER_ID:MSVC>:/W4 /O2 /arch:AVX2 /fp:fast>
    )

    if (MONOFIT_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${NAME} PRIVATE -fsanitize=address,undefined)
        target_link_libraries(${NAME} PRIVATE -fsanitize=address,undefined)
    endif ()

    add_test(NAME ${NAME} COMMAND ${NAME})
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
    target_compile_features(${NAME} PRIVATE cxx_std_20)

    target_compile_options(${NAME} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native -ffast-math>
            $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2>
    )
    # not a CTest test by default
endfunction()

add_polyfit_test(test_horner test_horner.cpp)
add_polyfit_benchmark(bench_horner bench_horner.cpp)

add_polyfit_test(test1D test_1D.cpp)
add_polyfit_benchmark(bench1D bench_1D.cpp)


add_polyfit_test(testMany test_many.cpp)
add_polyfit_benchmark(benchMany bench_many.cpp)

add_polyfit_test(test_ND test_ND.cpp)
add_polyfit_benchmark(bench_ND bench_ND.cpp)

add_polyfit_test(test_general test_general.cpp)
enable_testing()
