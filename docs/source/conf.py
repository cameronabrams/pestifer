# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'pestifer'
copyright = '2023-2025, Cameron F. Abrams'
author = 'cfa22@drexel.edu'

release = '1.15'
version = '1.15.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinxcontrib.mermaid',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

# html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
# html_theme_options = {
#     "footer_icons": [
#         {
#             "name": "Test",
#             "url": "https://example.com",
#             "html": """
#                 <svg width="24" height="24" xmlns="http://www.w3.org/2000/svg">
#                     <rect width="24" height="24" style="fill:blue;"/>
#                 </svg>
#             """,
#             "class": "example-icon",
#         }
#     ]
# }

html_theme_options = {
    "light_css_variables": {
        "color-icon": "#000000"  # Black for light mode
    },
    "dark_css_variables": {
        "color-icon": "#FFFFFF"  # White for dark mode
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/cameronabrams/pestifer",
            "html": """
    <svg role="img" width="24" height="24" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg" fill="currentColor">
        <title>GitHub</title>
        <path d="M12 0C5.37 0 0 5.37 0 12c0 5.3 3.44 9.8 8.21 11.39.6.11.82-.26.82-.58v-2.03c-3.34.73-4.04-1.61-4.04-1.61-.54-1.38-1.33-1.75-1.33-1.75-1.09-.75.08-.74.08-.74 1.2.08 1.83 1.23 1.83 1.23 1.07 1.83 2.81 1.3 3.5.99.11-.77.42-1.3.76-1.6-2.67-.3-5.47-1.34-5.47-5.98 0-1.32.47-2.4 1.24-3.24-.12-.3-.54-1.51.12-3.14 0 0 1.01-.32 3.3 1.23a11.38 11.38 0 0 1 3 0c2.28-1.55 3.3-1.23 3.3-1.23.66 1.63.24 2.84.12 3.14.77.84 1.24 1.92 1.24 3.24 0 4.65-2.8 5.68-5.47 5.98.43.37.81 1.1.81 2.22v3.29c0 .32.22.69.82.58C20.56 21.8 24 17.3 24 12c0-6.63-5.37-12-12-12z"/>
    </svg>
""",

            "class": "github-icon",
        },
        {
            "name": "LinkedIn",
            "url": "https://linkedin.com/in/cameron-abrams-b0143398",
            "html": """
                <svg role="img" width="24" height="24" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg" fill="currentColor">
                    <title>LinkedIn</title>
                    <path d="M20.447 20.452h-3.554v-5.403c0-1.288-.025-2.945-1.796-2.945-1.796 0-2.071 1.4-2.071 2.847v5.501h-3.554V9.001h3.414v1.561h.05c.475-.9 1.637-1.797 3.368-1.797 3.599 0 4.262 2.368 4.262 5.446v6.241zM5.337 7.433c-1.144 0-2.072-.93-2.072-2.072 0-1.142.928-2.07 2.072-2.07 1.142 0 2.07.928 2.07 2.07 0 1.144-.928 2.072-2.07 2.072zM6.814 20.452H3.859V9.001h2.955v11.451zM22.225 0H1.771C.792 0 0 .774 0 1.729v20.543C0 23.225.792 24 1.771 24h20.451C23.208 24 24 23.225 24 22.272V1.729C24 .774 23.208 0 22.225 0z"/>
                </svg>
            """,
            "class": "linkedin-icon",
        }
    ],
}

html_static_path = ['_static']

# -- Options for EPUB output
epub_show_urls = 'footnote'

mermaid_params = ['--theme', 'dark', '--width', '600']

def setup(app):
    app.add_css_file("css/custom.css")
    