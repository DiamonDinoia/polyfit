option(POLYFIT_GENERATE_DOCS "Generate Doxygen + Sphinx documentation" OFF)

if (NOT POLYFIT_GENERATE_DOCS)
    return()
endif ()

find_package(Doxygen REQUIRED)
find_program(SPHINX_BUILD sphinx-build REQUIRED)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile"
    @ONLY
)

add_custom_target(doxygen
    COMMAND Doxygen::doxygen "${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating Doxygen XML"
    VERBATIM
)

add_custom_target(sphinx
    COMMAND ${CMAKE_COMMAND} -E env
        DOXYGEN_XML_OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/docs/xml"
        ${SPHINX_BUILD} -b html
        "${CMAKE_CURRENT_SOURCE_DIR}/docs"
        "${CMAKE_CURRENT_BINARY_DIR}/docs/html"
    DEPENDS doxygen
    COMMENT "Generating Sphinx HTML docs"
    VERBATIM
)

add_custom_target(docs DEPENDS sphinx
    COMMENT "Build complete: open build/docs/html/index.html"
)
