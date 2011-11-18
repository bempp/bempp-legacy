# Additional checks needed to build the module
AC_DEFUN([DUNE_FOAMGRID_CHECKS])

# Additional checks needed to find the module
AC_DEFUN([DUNE_FOAMGRID_CHECK_MODULE]),[
  DUNE_CHECK_MODULES([dune-foamgrid], [dune-foamgrid/foamgrid.hh])
])
