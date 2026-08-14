#pragma once

// Single definition site for `stdex`. CMake sets POLYFIT_USE_STD_MDSPAN so every
// TU in a build agrees; the __cpp_lib_mdspan probe below is the fallback for use
// without CMake.

#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif

#if !defined(POLYFIT_USE_STD_MDSPAN)
#  if defined(__cpp_lib_mdspan) && __cpp_lib_mdspan >= 202207L
#    define POLYFIT_USE_STD_MDSPAN 1
#  else
#    define POLYFIT_USE_STD_MDSPAN 0
#  endif
#endif

#if POLYFIT_USE_STD_MDSPAN
#include <mdspan>
namespace stdex = std;
#else
#include <experimental/mdspan>
namespace stdex = std::experimental;
#endif
