# Configuration file for the Sphinx documentation builder.
#
# -- Project information -----------------------------------------------------

project = 'LFRic Core'
author = 'Core Capability Development Team'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'vcs_pageview_mode': '',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom CSS override file.
html_style = 'css/custom.css'

# Provides the Edit on GitHub link in the generated docs.
html_context = {
    "display_github": True,
    "github_user": "MetOffice",
    "github_repo": "lfric_core",
    "github_version": "main",
    "conf_py_path": "/documentation/source/"
}

# Enable numbered references to e.g. figures.
#
numfig = True

