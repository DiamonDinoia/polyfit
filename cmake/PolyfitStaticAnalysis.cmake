include_guard(GLOBAL)

option(POLYFIT_ENABLE_CLANG_TIDY "Enable clang-tidy diagnostics" OFF)
option(POLYFIT_CLANG_TIDY_CHECKS "Override clang-tidy checks" "")
option(POLYFIT_CLANG_TIDY_WARNINGS_AS_ERRORS "Treat clang-tidy warnings as errors" ON)
option(POLYFIT_ENABLE_CPPCHECK "Enable cppcheck diagnostics" OFF)
option(POLYFIT_CPPCHECK_OPTIONS "Additional cppcheck options" "--enable=warning,style,performance,portability")

function(polyfit_configure_static_analysis target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "polyfit_configure_static_analysis called with non-existent target '${target}'")
  endif()

  # Skip analysis for test and benchmark targets — only lint examples and library headers.
  get_target_property(_src_dir ${target} SOURCE_DIR)
  string(FIND "${_src_dir}" "${CMAKE_SOURCE_DIR}/tests" _is_test)
  if(NOT _is_test EQUAL -1)
    return()
  endif()

  if(POLYFIT_ENABLE_CLANG_TIDY)
    find_program(_clang_tidy_exe NAMES clang-tidy clang-tidy-21 clang-tidy-20 clang-tidy-19 clang-tidy-18 clang-tidy-17 clang-tidy-16)
    if(_clang_tidy_exe)
      set(_clang_tidy_cmd "${_clang_tidy_exe}")
      set(_clang_tidy_checks "${POLYFIT_CLANG_TIDY_CHECKS}")
      if(NOT _clang_tidy_checks)
        set(_clang_tidy_checks
            "-*,clang-analyzer-*,bugprone-*,-bugprone-easily-swappable-parameters,-bugprone-exception-escape,performance-*,misc-const-correctness,misc-misplaced-const,misc-redundant-expression,misc-unused-parameters,readability-inconsistent-declaration-parameter-name,readability-misleading-indentation,readability-redundant-control-flow,readability-redundant-member-init,readability-redundant-preprocessor,readability-simplify-boolean-expr,readability-static-accessed-through-instance,modernize-use-nodiscard")
      endif()
      list(APPEND _clang_tidy_cmd "-checks=${_clang_tidy_checks}")
      list(APPEND _clang_tidy_cmd "-header-filter=^${CMAKE_SOURCE_DIR}/include/polyfit")
      list(APPEND _clang_tidy_cmd "--extra-arg=-fsyntax-only")
      # clang-tidy replays the compile command with Clang, so GCC-only warnings
      # from the normal warning profile must not become analysis failures.
      list(APPEND _clang_tidy_cmd "--extra-arg=-Wno-unknown-warning-option")
      if(POLYFIT_CLANG_TIDY_WARNINGS_AS_ERRORS)
        list(APPEND _clang_tidy_cmd "-warnings-as-errors=*")
      endif()
      set_property(TARGET ${target} PROPERTY CXX_CLANG_TIDY "${_clang_tidy_cmd}")
    else()
      message(WARNING "POLYFIT_ENABLE_CLANG_TIDY is ON but clang-tidy was not found on PATH")
    endif()
  endif()

  if(POLYFIT_ENABLE_CPPCHECK)
    find_program(_cppcheck_exe NAMES cppcheck)
    if(_cppcheck_exe)
      # -D__cppcheck__: force the guard active so PF_IF_CONSTEVAL uses the
      # C++20-compatible form (cppcheck's parser cannot handle `if consteval`).
      set(_cppcheck_cmd "${_cppcheck_exe}" "--inline-suppr" "-D__cppcheck__")
      if(POLYFIT_CPPCHECK_OPTIONS)
        list(APPEND _cppcheck_cmd "${POLYFIT_CPPCHECK_OPTIONS}")
      endif()
      set_property(TARGET ${target} PROPERTY CXX_CPPCHECK "${_cppcheck_cmd}")
    else()
      message(WARNING "POLYFIT_ENABLE_CPPCHECK is ON but cppcheck was not found on PATH")
    endif()
  endif()
endfunction()
