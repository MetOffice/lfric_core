# Configuration file for the Sphinx documentation builder.
#
# -- imports -----------------------------------------------------------------
import os
# -- Project information -----------------------------------------------------

project = 'LFRic Core'
author = 'Core Capability Development Team'
copyright = '2024 Met Office. All rights reserved'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_sitemap',
    'sphinx_design',
    'sphinx.ext.intersphinx'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'
# html_title = "LFRic Core"

# Generate the sitemap info, this will need updating when we have versioned docs
html_baseurl = os.environ.get("SPHINX_HTML_BASE_URL", "https://metoffice.github.io/lfric_core")
sitemap_locales = [None]
sitemap_url_scheme = "{link}"

# Hide the link which shows the rst markup
html_show_sourcelink = False

html_theme_options = {
    "navigation_with_keys": True,
    "use_edit_page_button": True,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_align": "content",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/MetOffice/lfric_core",
            "icon": "fa-brands fa-github"
        },
        {
            "name": "GitHub Discussions",
            "url": "https://github.com/MetOffice/simulation-systems/discussions",
            "icon": "far fa-comments",
        }
    ],
    "logo": {
        "text": "LFRic Core",
        "image_light": "_static/MO_SQUARE_black_mono_for_light_backg_RBG.png",
        "image_dark": "_static/MO_SQUARE_for_dark_backg_RBG.png",
    },
    "secondary_sidebar_items": {
        "**/*": ["page-toc", "edit-this-page", "show-glossary", "show-faqs"],
        "index": [],
    },
    "footer_start": ["crown-copyright"],
    "footer_center": ["sphinx-version"],
    "footer_end": ["theme-version"],
    "primary_sidebar_end": []
}

html_sidebars = {
    "index": []
}

# Provides the Edit on GitHub link in the generated docs.
html_context = {
    "display_github": True,
    "github_user": "MetOffice",
    "github_repo": "lfric_core",
    "github_version": "main",
    "doc_path": "documentation/source"
}

# Enable numbered references to e.g. figures.
#
numfig = True

# Options for intersphinx mapping extension
# To discover the available objects in the target project (e.g. psyclone), run:
# python -m sphinx.ext.intersphinx https://psyclone.readthedocs.io/en/stable/objects.inv

intersphinx_mapping = {
    'psyclone': ('https://psyclone.readthedocs.io/en/stable/', None),
    'simsys': ('https://metoffice.github.io/simulation-systems/', None)
}

