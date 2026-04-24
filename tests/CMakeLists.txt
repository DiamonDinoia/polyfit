cmake_minimum_required(VERSION 3.16)
project(polyfit_tests LANGUAGES C CXX)

# -----------------------------------------------------------------------------
# Basic config
# -----------------------------------------------------------------------------
set(CMAKE_CXX_EXTENSIONS OFF)

include(PolyfitWarnings)
include(PolyfitSanitizers)
include(PolyfitStaticAnalysis)

# -----------------------------------------------------------------------------
# Fetch GoogleTest & Benchmark
# -----------------------------------------------------------------------------
CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION 1.17.0
        OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt ON"
        "BUILD_GMOCK OFF"
)


CPMAddPackage(
        nanobench
        GITHUB_REPOSITORY martinus/nanobench
        GIT_TAG v4.3.11
)

option(POLYFIT_ENABLE_BENCH_FAST_MATH "Enable fast-math for benchmark targets" OFF)
option(POLYFIT_ENABLE_BENCH_NATIVE_TUNING "Enable -march=native for benchmark targets" ${POLYFIT_ENABLE_NATIVE_TUNING})
option(POLYFIT_BUILD_BENCHMARKS "Build benchmark executables" ON)

function(polyfit_apply_dev_profile target)
    if (POLYFIT_STRICT_WARNINGS)
        polyfit_enable_warnings(${target})
    endif ()

    polyfit_enable_sanitizers(${target})
    polyfit_configure_static_analysis(${target})

    if (POLYFIT_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif ()
endfunction()

function(polyfit_apply_bench_profile target)
    if (POLYFIT_STRICT_WARNINGS)
        polyfit_enable_warnings(${target})
    endif ()
endfunction()

function(polyfit_configure_benchmark_target target std)
    target_link_libraries(${target}
            PRIVATE
            polyfit::polyfit
            nanobench::nanobench
    )
    set_property(TARGET ${target} PROPERTY CXX_STANDARD ${std})
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${target} PROPERTY CXX_EXTENSIONS OFF)
    target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-O3>
            $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<NOT:$<CONFIG:Debug>>>:/O2>
            $<$<AND:$<BOOL:${POLYFIT_ENABLE_BENCH_FAST_MATH}>,$<CXX_COMPILER_ID:GNU,Clang>>:-ffast-math>
            $<$<AND:$<BOOL:${POLYFIT_ENABLE_BENCH_NATIVE_TUNING}>,$<CXX_COMPILER_ID:GNU,Clang>>:-march=native>
            $<$<AND:$<BOOL:${POLYFIT_ENABLE_BENCH_FAST_MATH}>,$<CXX_COMPILER_ID:MSVC>>:/fp:fast>
    )
    polyfit_apply_bench_profile(${target})
endfunction()

# -----------------------------------------------------------------------------
# Helpers to add unit-test targets for multiple C++ standards
# -----------------------------------------------------------------------------
# Standards to test/build: by default start at C++17 and add every later
# standard the compiler supports. CI jobs may override this cache variable to
# trim expensive configurations, for example Valgrind-only memory checks.
set(POLYFIT_TEST_STANDARDS "" CACHE STRING
        "Semicolon-separated C++ standards to build for tests and benchmarks; empty means auto-detect from C++17 upward")
if (POLYFIT_TEST_STANDARDS STREQUAL "")
    set(POLYFIT_TEST_STANDARDS 17)
    include(CheckCXXCompilerFlag)
    foreach(_std IN ITEMS 20 23 26)
        check_cxx_compiler_flag("-std=c++${_std}" _POLYFIT_CXX${_std})
        if (_POLYFIT_CXX${_std})
            list(APPEND POLYFIT_TEST_STANDARDS ${_std})
        endif ()
    endforeach()
endif ()
message(STATUS "Polyfit test standards: ${POLYFIT_TEST_STANDARDS}")

function(add_polyfit_test_for_std BASENAME STD)
    set(TGT "${BASENAME}_cxx${STD}")
    add_executable(${TGT} ${ARGN})
    target_link_libraries(${TGT}
            PRIVATE polyfit::polyfit gtest_main
    )
    # Per-target standard selection
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD ${STD})
    set_property(TARGET ${TGT} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${TGT} PROPERTY CXX_EXTENSIONS OFF)
    polyfit_apply_dev_profile(${TGT})

    add_test(NAME ${TGT} COMMAND ${TGT})
    set_property(GLOBAL APPEND PROPERTY POLYFIT_TEST_TARGETS ${TGT})
endfunction()

function(add_polyfit_test_matrix NAME)
    foreach(STD IN LISTS POLYFIT_TEST_STANDARDS)
        add_polyfit_test_for_std(${NAME} ${STD} ${ARGN})
    endforeach()
endfunction()

# Benchmarks matrix across standards (where supported)
function(add_polyfit_benchmark_for_std BASENAME STD)
    set(TGT "${BASENAME}_cxx${STD}")
    add_executable(${TGT} ${ARGN})
    polyfit_configure_benchmark_target(${TGT} ${STD})
    set_property(GLOBAL APPEND PROPERTY POLYFIT_BENCHMARK_TARGETS_${STD} ${TGT})
endfunction()

function(add_polyfit_benchmark_matrix NAME)
    foreach(STD IN LISTS POLYFIT_TEST_STANDARDS)
        add_polyfit_benchmark_for_std(${NAME} ${STD} ${ARGN})
    endforeach()
endfunction()

add_polyfit_test_matrix(test_horner test_horner.cpp)
if (POLYFIT_BUILD_BENCHMARKS)
    add_polyfit_benchmark_matrix(bench_horner bench_horner.cpp)
    add_polyfit_benchmark_matrix(bench_unroll_sweep bench_unroll_sweep.cpp)
endif ()

add_polyfit_test_matrix(test1D_basic      test_1D_basic.cpp)
add_polyfit_test_matrix(test1D_truncation test_1D_truncation.cpp)
add_polyfit_test_matrix(test1D_accuracy   test_1D_accuracy.cpp)
add_polyfit_test_matrix(test1D_fusion     test_1D_fusion.cpp)
if (POLYFIT_BUILD_BENCHMARKS)
    add_polyfit_benchmark_matrix(bench1D bench_1D.cpp)
endif ()


add_polyfit_test_matrix(testMany test_many.cpp)
if (POLYFIT_BUILD_BENCHMARKS)
    add_polyfit_benchmark_matrix(benchMany bench_many.cpp)
endif ()

add_polyfit_test_matrix(test_ND test_ND_accuracy.cpp test_ND_api.cpp test_ND_fusion.cpp)
if (POLYFIT_BUILD_BENCHMARKS)
    add_polyfit_benchmark_matrix(bench_ND bench_ND.cpp)
endif ()

add_polyfit_test_matrix(test_general test_general.cpp)
add_polyfit_test_matrix(test_copy_move test_copy_move.cpp)
add_polyfit_test_matrix(test_scaling_alignment test_scaling_alignment.cpp)
add_polyfit_test_matrix(test_utils_bjorck_newton test_utils_bjorck_newton.cpp)
add_polyfit_test_matrix(test_utils_linspace test_utils_linspace.cpp)
add_polyfit_test_matrix(test_func_eval_many_padding test_func_eval_many_padding.cpp)
add_polyfit_test_matrix(test_edge_cases test_edge_cases.cpp)
add_polyfit_test_matrix(test_constexpr test_constexpr.cpp)
add_polyfit_test_matrix(test_feature_macros test_feature_macros.cpp)
add_polyfit_test_matrix(test_portable_trig test_portable_trig.cpp)

if (POLYFIT_BUILD_BENCHMARKS)
    foreach(STD IN LISTS POLYFIT_TEST_STANDARDS)
        get_property(_polyfit_bench_targets GLOBAL PROPERTY POLYFIT_BENCHMARK_TARGETS_${STD})
        if (_polyfit_bench_targets)
            set(_runner "run_benchmarks_cxx${STD}")
            set(_runner_commands
                    COMMAND ${CMAKE_COMMAND} -E echo "Benchmark note: use a tuned system for stable numbers."
                    COMMAND ${CMAKE_COMMAND} -E echo "Recommended: pyperf system tune, performance governor, and native tuning if desired."
            )
            foreach(_bench_target IN LISTS _polyfit_bench_targets)
                list(APPEND _runner_commands COMMAND $<TARGET_FILE:${_bench_target}>)
            endforeach ()
            add_custom_target(${_runner}
                    ${_runner_commands}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Running polyfit benchmarks for C++${STD}"
                    VERBATIM
            )
            add_dependencies(${_runner} ${_polyfit_bench_targets})
            if (STD EQUAL 20)
                add_custom_target(run_benchmarks DEPENDS ${_runner})
            endif ()
        endif ()
    endforeach ()
endif ()

if (POLYFIT_ENABLE_COVERAGE)
    find_program(LCOV_EXECUTABLE lcov)
    find_program(GENHTML_EXECUTABLE genhtml)
    find_program(GCOVR_EXECUTABLE gcovr)

    get_property(_polyfit_all_tests GLOBAL PROPERTY POLYFIT_TEST_TARGETS)

    if (LCOV_EXECUTABLE AND GENHTML_EXECUTABLE)
        set(LCOV_INFO "${CMAKE_BINARY_DIR}/coverage.info")
        set(LCOV_FILTERED "${CMAKE_BINARY_DIR}/coverage.filtered.info")
        set(COVERAGE_DIR "${CMAKE_BINARY_DIR}/coverage")

        add_custom_target(coverage
                COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR} --output-on-failure
                COMMAND ${LCOV_EXECUTABLE} --capture --directory ${CMAKE_BINARY_DIR} --output-file ${LCOV_INFO} --ignore-errors inconsistent,unused,mismatch
                COMMAND ${LCOV_EXECUTABLE} --remove ${LCOV_INFO} "/usr/*" "*/_deps/*" --output-file ${LCOV_FILTERED} --ignore-errors inconsistent,unused,mismatch
                COMMAND ${GENHTML_EXECUTABLE} --ignore-errors mismatch -o ${COVERAGE_DIR} ${LCOV_FILTERED}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                COMMENT "Running tests and generating coverage report (lcov+genhtml)"
                VERBATIM
        )
    elseif (GCOVR_EXECUTABLE)
        add_custom_target(coverage
                COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR} --output-on-failure
                COMMAND ${GCOVR_EXECUTABLE} -r ${CMAKE_SOURCE_DIR} --exclude ".*/_deps/.*" --exclude "/usr/.*" --gcov-ignore-errors=no_working_dir_found --html --html-details -o ${CMAKE_BINARY_DIR}/coverage-report.html
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                COMMENT "Running tests and generating coverage report (gcovr)"
                VERBATIM
        )
    else ()
        add_custom_target(coverage
                COMMAND ${CMAKE_COMMAND} -E echo "Coverage tools missing. Install 'lcov'+'genhtml' or 'gcovr'."
        )
    endif ()

    if (TARGET coverage AND _polyfit_all_tests)
        add_dependencies(coverage ${_polyfit_all_tests})
    endif ()
endif ()

enable_testing()
