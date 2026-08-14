#pragma once

// Single definition site for `stdex`. <version> first, so the alias never
// depends on include order. The threshold matches POLYFIT_HAS_NATIVE_MDSPAN
// in CMakeLists.txt; both decide the kokkos/mdspan package identically.

#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif

#if defined(__cpp_lib_mdspan) && __cpp_lib_mdspan >= 202207L
#include <mdspan>
namespace stdex = std;
#else
#include <experimental/mdspan>
namespace stdex = std::experimental;
#endif
