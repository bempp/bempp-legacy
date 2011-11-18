
SUBDIRS = doxygen

CURDIR = doc

BASEDIR = ..

docdir=$(datadir)/doc/dune-foamgrid

include $(top_srcdir)/am/webstuff

CLEANFILES = $(PAGES)

if ! BUILD_DOCS
# add tag to notify that dist has been build without documentation
dist-hook:
        echo "# No documentation included in distribution! " > $(distdir)/$(DOCUMENTATION_TAG_FILE)
endif

include $(top_srcdir)/am/global-rules
