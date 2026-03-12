include_guard(GLOBAL)

option(POLYFIT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)

function(polyfit_enable_warnings target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "polyfit_enable_warnings called with non-existent target '${target}'")
  endif()

  get_target_property(_target_type "${target}" TYPE)
  if(_target_type STREQUAL "INTERFACE_LIBRARY")
    set(_scope INTERFACE)
  else()
    set(_scope PRIVATE)
  endif()

  set(_clang_like $<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>)
  set(_gnu $<CXX_COMPILER_ID:GNU>)
  set(_gnu_or_clang $<OR:${_gnu},${_clang_like}>)
  set(_msvc $<CXX_COMPILER_ID:MSVC>)
  set(_lang_is_cxx $<COMPILE_LANGUAGE:CXX>)

  # Warnings enabled for both GCC and Clang.
  set(_warnings_common
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wdouble-promotion
    -Wold-style-cast
    -Wnon-virtual-dtor
    -Woverloaded-virtual
    -Wcast-align
    -Wunused
    -Wimplicit-fallthrough
    -Wformat=2
  )

  # Clang-only warnings.
  # - Wnull-dereference: GCC's IPA produces false positives when deeply inlining
  #   standard library vector operations through template chains.
  # - Wno-gnu-zero-variadic-macro-arguments: GTest TYPED_TEST_SUITE uses GNU
  #   extension zero-arg variadic macros; suppress the pedantic warning.
  # - Wno-c++20-attribute-extensions: allow using standard attributes such as
  #   [[likely]] and [[no_unique_address]] as extensions in C++17 mode.
  set(_warnings_clang_only
    -Wnull-dereference
    -Wno-gnu-zero-variadic-macro-arguments
    -Wno-c++20-attribute-extensions
  )

  set(_warnings_msvc
    /W4
    /permissive-
    /Zc:__cplusplus  # report correct __cplusplus value instead of always 199711L
    /external:anglebrackets  # treat angle-bracket includes as external
    /external:W0             # suppress all warnings in external headers
    /wd4702  # unreachable code — false positive with if constexpr branches
    /wd4723  # potential divide by 0 — false positive in interval scaling
    /w14242
    /w14254
    /w14263
    /w14265
    /w14287
    /we4289
    /w14296
    /w14311
    /w14545
    /w14546
    /w14547
    /w14549
    /w14555
    /w14619
    /w14640
    /w14826
    /w14905
    /w14906
    /w14928
  )

  # Additional GCC-only warnings (checked at configure time; applied only to GCC).
  set(_additional_warnings
    -Wduplicated-cond
    -Wlogical-op
    -Wuseless-cast
    -Winit-self
    -Wmissing-include-dirs
    -Wredundant-decls
    -Wmisleading-indentation
    -Wsuggest-override
  )

  set(_compile_options)
  foreach(_flag IN LISTS _warnings_common)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_gnu_or_clang}>:${_flag}>)
  endforeach()
  foreach(_flag IN LISTS _warnings_clang_only)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_clang_like}>:${_flag}>)
  endforeach()
  foreach(_flag IN LISTS _warnings_msvc)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_msvc}>:${_flag}>)
  endforeach()

  include(CheckCXXCompilerFlag)
  foreach(_flag IN LISTS _additional_warnings)
    string(MAKE_C_IDENTIFIER "POLYFIT_HAS${_flag}" _flag_id)
    check_cxx_compiler_flag("${_flag}" "${_flag_id}")
    if(${_flag_id})
      # Apply additional (GCC-specific) warnings only to GCC, not Clang/MSVC,
      # to avoid "unknown warning option" errors in clang-tidy and Clang builds.
      list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_gnu}>:${_flag}>)
    endif()
  endforeach()

  if(POLYFIT_WARNINGS_AS_ERRORS)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_gnu_or_clang}>:-Werror>)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_msvc}>:/WX>)
  endif()

  target_compile_options(${target} ${_scope} ${_compile_options})
endfunction()
