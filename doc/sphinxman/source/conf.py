# -*- coding: utf-8 -*-
#
# Psithon documentation build configuration file, created by
# sphinx-quickstart on Sun Feb 12 04:25:25 2012.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('../../../lib/python/'))
sys.path.insert(0, os.path.abspath('../../../lib/databases/'))
sys.path.insert(0, os.path.abspath('../../../tests/'))
sys.path.insert(0, os.path.abspath('.'))

import psi4_sptheme as psp

# very bad practice- hardwiring object dir so can get at gitversion.h
objdir = 'objdir'

# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.1'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
# 'psi4_sptheme.ext.autodoc_sections',
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.pngmath', 'sphinx.ext.viewcode', 
    'sphinx.ext.extlinks', 'psi4_sptheme.ext.index_styling', 'psi4_sptheme.ext.psidomain']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'PSI4'
copyright = u'2012, Psi4 Project'

def get_gitversion():
    file = '../../../' + objdir + '/src/bin/psi4/gitversion.h'

    try: 
        f = open(file)
    except IOError:
        return ''
    vers = f.readline().split()[2]
    f.close() 
    print 'GIT_VERSION', vers[1:7]
    return vers[1:7]

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = 'devel'
# The full version, including alpha/beta/rc tags.
#release = version + '.' + get_gitversion()
release = get_gitversion()

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
#exclude_patterns = []
exclude_patterns = ['template_index.rst',          'template_appendices.rst',
                    'autodoc_index_html.rst',      'autodoc_appendices_html.rst', 
                    'autodoc_index_htmluser.rst',  'autodoc_appendices_htmluser.rst', 
                    'autodoc_index_htmlprog.rst',  'autodoc_appendices_htmlprog.rst', 
                    'autodoc_index_latexuser.rst', 'autodoc_appendices_latexuser.rst', 
                    'autodoc_index_latexprog.rst', 'autodoc_appendices_latexprog.rst',
                    'abbr_accents.rst',
                    'autodoc_abbr_options_c.rst',
                    'autodoc_abbr_options_plugins.rst']


# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False
show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'psi4'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {"roottarget": "index" }

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [psp.get_theme_dir()]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None
html_title = 'PSI4 [' + release + '] documentation'

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = "PSI4_3.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {'**': ['localtoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'Psi4doc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'Psi4.tex', u'Psi4 Documentation',
   u'Psi4 Project', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'psithon', u'Psithon Documentation',
     [u'Psi4 Project'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'Psithon', u'Psithon Documentation',
   u'Psi4 Project', 'Psithon', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# Abbreviations
rst_epilog = """
.. |PSIfour| replace:: PSI4
.. |psirc| replace:: ``~/.psi4rc``
.. |dl| replace:: :math:`\Rightarrow`
.. |dr| replace:: :math:`\Leftarrow`
.. |kcalpermol| replace:: kcal mol\ :sup:`-1`
.. |Angstrom| replace:: |AA|\ ngstr\ |o_dots|\ m
.. |--| unicode:: U+02013 .. en dash
   :trim:
.. |---| unicode:: U+02014 .. em dash
   :trim:
.. |w---w| unicode:: U+02014 .. em dash
.. include:: /abbr_accents.rst
"""

# Logo at top of page
rst_prolog = """
.. image:: /psi4banner.png
   :width: 100 %
   :alt: PSI4 Project Logo
"""

# This option, from sphinx.ext.extlinks, allows linking to source on trac with :source:`lib/python/driver/py`
extlinks = {'source': ('http://sirius.chem.vt.edu/trac/browser/%s', 'psi4/'),
            'srcsample': ('http://sirius.chem.vt.edu/trac/browser/samples/%s/input.dat', ''),
            'srcbasis': ('http://sirius.chem.vt.edu/trac/browser/lib/basis/%s.gbs', '') }

# corrects baseline for inline math + other custom stuff for inline math, such as non-default math fonts etc.
pngmath_latex_preamble=r'\usepackage[active]{preview}'
pngmath_use_preview=True



