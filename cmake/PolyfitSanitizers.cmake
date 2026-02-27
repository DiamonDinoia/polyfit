include_guard(GLOBAL)

option(POLYFIT_ENABLE_SANITIZERS "Master switch to enable sanitizers" OFF)
option(POLYFIT_ENABLE_ASAN "Enable AddressSanitizer" ${POLYFIT_ENABLE_SANITIZERS})
option(POLYFIT_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" ${POLYFIT_ENABLE_SANITIZERS})

function(polyfit_enable_sanitizers target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "polyfit_enable_sanitizers called with non-existent target '${target}'")
  endif()

  if(NOT (POLYFIT_ENABLE_ASAN OR POLYFIT_ENABLE_UBSAN))
    return()
  endif()

  if(MSVC)
    message(WARNING "POLYFIT sanitizers are currently only configured for GCC/Clang-compatible toolchains")
    return()
  endif()

  get_target_property(_target_type "${target}" TYPE)
  if(_target_type STREQUAL "INTERFACE_LIBRARY")
    set(_scope INTERFACE)
  else()
    set(_scope PRIVATE)
  endif()

  set(_sanitize_flags)
  if(POLYFIT_ENABLE_ASAN)
    list(APPEND _sanitize_flags -fsanitize=address -fno-omit-frame-pointer)
  endif()
  if(POLYFIT_ENABLE_UBSAN)
    list(APPEND _sanitize_flags -fsanitize=undefined)
  endif()
  list(REMOVE_DUPLICATES _sanitize_flags)

  set(_gnu_or_clang $<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>)
  set(_lang_is_cxx $<COMPILE_LANGUAGE:CXX>)

  set(_compile_options)
  set(_link_options)
  foreach(_flag IN LISTS _sanitize_flags)
    list(APPEND _compile_options $<$<AND:${_lang_is_cxx},${_gnu_or_clang}>:${_flag}>)
    list(APPEND _link_options $<$<${_gnu_or_clang}>:${_flag}>)
  endforeach()

  target_compile_options(${target} ${_scope} ${_compile_options})
  target_link_options(${target} ${_scope} ${_link_options})
endfunction()
