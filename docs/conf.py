import os
project = "polyfit"
copyright = "2025, polyfit contributors"
author = "polyfit contributors"
release = "1.0.0"

extensions = [
    "breathe",
    "exhale",
    "myst_parser",
    "sphinx_rtd_theme",
]

html_theme = "sphinx_rtd_theme"

# Breathe: point at Doxygen XML output
_doxygen_xml = os.environ.get(
    "DOXYGEN_XML_OUTPUT",
    os.path.join(os.path.dirname(__file__), "../build/docs/xml"),
)
breathe_projects = {"polyfit": _doxygen_xml}
breathe_default_project = "polyfit"

# Exhale: auto-generate API tree from Breathe
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "API Reference",
    "doxygenStripFromPath": os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "include", "polyfit")),
    "createTreeView": True,
    "exhaleExecutesDoxygen": False,
    # Exclude internal implementation detail namespace
    "unabridgedOrphanKinds": {"namespace"},
}

# MyST: parse .md files
myst_enable_extensions = ["colon_fence", "deflist"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "API.md",
    "api/function_polyeval_8hpp_1a15ae258635ac1bb6ac5a10da6d03bf09.rst",
    "api/function_polyeval_8hpp_1a1b97c0ec59989286dda3336ec1a3dc38.rst",
    "api/function_polyeval_8hpp_1a4c21a6ec3ce2229ddfbd84898b9520be.rst",
    "api/function_polyeval_8hpp_1ab79217b318d2b17aae26986ea328ec12.rst",
    "api/function_polyeval_8hpp_1ac9707fe8b6105bc85135351dfd24bf7a.rst",
    "api/function_polyeval_8hpp_1ae71c154976d4f3007d0d78ba81d08e76.rst",
    "api/function_polyeval_8hpp_1aeabd8326d5cc248fbd6a8ae8d1098f40.rst",
    "api/function_polyeval_8hpp_1afa1c2905c64c53c7d718c72a49eec305.rst",
    "api/namespace_std.rst",
]

html_static_path = []
html_show_sourcelink = False
