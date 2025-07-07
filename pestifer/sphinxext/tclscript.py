"""
This Sphinx extension provides a directive to include Tcl scripts in documentation.
Written by ChatGPT, based on the original work by Cameron F. Abrams.
"""

import os
from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import ViewList
from sphinx.util.nodes import nested_parse_with_titles

class TclScriptDirective(Directive):
    required_arguments = 1
    has_content = False

    def run(self):
        script_path = os.path.normpath(self.arguments[0])

        # Repo root: two levels up from this file (docs/source/)
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        full_path = os.path.normpath(os.path.join(repo_root, script_path))

        if not os.path.isfile(full_path):
            return [nodes.paragraph(text=f"ERROR: File not found: {full_path}")]

        # Track this file as a dependency for rebuilds
        env = self.state.document.settings.env
        env.note_dependency(full_path)

        # Extract header comment block
        with open(full_path, 'r') as f:
            lines = f.readlines()

        header_lines = []
        for line in lines:
            if line.strip().startswith('##'):
                header_lines.append(line.lstrip('#').strip())
            elif line.strip() == '' or line.strip().startswith('#'):
                header_lines.append('')
            else:
                break

        # GitHub source link
        github_relpath = os.path.relpath(full_path, repo_root).replace(os.sep, '/')
        github_url = f"https://github.com/cameronabrams/pestifer/blob/main/{github_relpath}"

        # Create a top-level section
        section = nodes.section(ids=[os.path.basename(script_path)])
        section += nodes.title(text=os.path.basename(script_path))

        # Parse the header as reStructuredText
        if header_lines:
            viewlist = ViewList()
            for i, line in enumerate(header_lines):
                viewlist.append(line, source=full_path, offset=i)

            nested_parse_with_titles(self.state, viewlist, section)

        # Add source link as a paragraph (no nested title!)
        link_para = nodes.paragraph()
        link_para += nodes.reference(text='[source]', refuri=github_url)
        section += link_para

        return [section]
