
TESTPROGS = foamgrid-test

# which tests to run
TESTS = $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = -Wall -g

# The directory where the example grid files used for testing are located.
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(DUNE_GRID_ROOT)/doc/grids/\"

# define the programs
foamgrid_test_SOURCES = foamgrid-test.cc

include $(top_srcdir)/am/global-rules
