# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import pysp2

# -- Project information -----------------------------------------------------

project = 'PySP2'
copyright = '2022, Robert Jackson'
author = 'Robert Jackson'

# The full version, including alpha/beta/rc tags

version = pysp2.__version__
release = pysp2.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_gallery.gen_gallery',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
    'sphinx_design',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

sphinx_gallery_conf = {
    'examples_dirs': '../../examples',
    'gallery_dirs': 'source/auto_examples'
}

# Generate the API documentation when building
autoclass_content = 'both'
autosummary_generate = True
autosummary_imported_members = True

# Otherwise, the Return parameter list looks different from the Parameter list
napoleon_use_rtype = False
napoleon_use_ivar = True
napoleon_include_init_with_doc = False
napoleon_use_param = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['default.css']

html_js_files = ['doc_shared.js']

html_sidebars = {
    'contributors_guide': ['searchbox.html', 'sidebar-nav-bs.html'],
    'developers_guide': ['searchbox.html', 'sidebar-nav-bs.html'],
    'users_guide': ['searchbox.html', 'sidebar-nav-bs.html'],
    'examples': ['searchbox.html', 'sidebar-nav-bs.html'],
    'notebook-gallery': ['searchbox.html', 'sidebar-nav-bs.html']}

