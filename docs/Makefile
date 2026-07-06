# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
APIDOC_OUT    = $(SOURCEDIR)/api
APIDOC_SRC    = ../pestifer

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help apidoc config-ref html clean Makefile

config-ref:
	cd $(SOURCEDIR) && yclept make-doc ../../pestifer/schema/base.yaml --root config_ref --footer-style raw-html

apidoc:
	rm -rf $(APIDOC_OUT)
	sphinx-apidoc -f -M -e --implicit-namespaces --tocfile API -o $(APIDOC_OUT) $(APIDOC_SRC) $(APIDOC_SRC)/resources

html: apidoc
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR) $(APIDOC_OUT)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
