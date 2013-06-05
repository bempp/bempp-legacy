
TESTPROGS =  foamgrid-test \
	global-refine-test \
	local-refine-test

# which tests to run
TESTS = $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = -Wall -g

# The directories where the example grid files used for testing are located.
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(DUNE_GRID_ROOT)/doc/grids/\"
AM_CPPFLAGS+=-DDUNE_FOAMGRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

# define the programs
foamgrid_test_SOURCES = foamgrid-test.cc

global_refine_test_SOURCES = global-refine-test.cc

local_refine_test_SOURCES = local-refine-test.cc

include $(top_srcdir)/am/global-rules
