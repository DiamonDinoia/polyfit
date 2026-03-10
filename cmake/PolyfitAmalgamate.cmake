find_program(Python3_EXECUTABLE NAMES python3 python)

if (NOT Python3_EXECUTABLE)
    message(STATUS "python3 not found — 'amalgamate' target unavailable")
    return()
endif ()

set(_AMALG_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/scripts/amalgamate.py")
set(_AMALG_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/include/polyfit/polyeval.hpp")
set(_AMALG_OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/include/polyfit/internal/polyfit_amalgamated.hpp")

# Glob the source headers so the target re-runs when any of them change.
file(GLOB_RECURSE _AMALG_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/polyfit/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/polyfit/*.hpp"
)
list(REMOVE_ITEM _AMALG_SOURCES "${_AMALG_OUTPUT}")

add_custom_command(
    OUTPUT "${_AMALG_OUTPUT}"
    COMMAND "${Python3_EXECUTABLE}" "${_AMALG_SCRIPT}"
        --root "${CMAKE_CURRENT_SOURCE_DIR}"
        --input "include/polyfit/polyeval.hpp"
        --output "include/polyfit/internal/polyfit_amalgamated.hpp"
    DEPENDS "${_AMALG_SCRIPT}" ${_AMALG_SOURCES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Generating single-header amalgamation: include/polyfit/internal/polyfit_amalgamated.hpp"
    VERBATIM
)

add_custom_target(amalgamate DEPENDS "${_AMALG_OUTPUT}"
    COMMENT "Single-header written to include/polyfit/internal/polyfit_amalgamated.hpp"
)
