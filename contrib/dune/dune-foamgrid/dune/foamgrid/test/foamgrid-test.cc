#include <config.h>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>


int main (int argc, char *argv[]) try
{
    // dimworld == 2
    //FoamGrid<2>* grid2d = make2DHybridTestGrid<FoamGrid<2> >();
    // paths to gmsh test files
    const std::string dune_grid_path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    std::auto_ptr<FoamGrid<2> > grid2d( GmshReader<FoamGrid<2> >::read( dune_grid_path + "curved2d.msh", false, false ) );

        
    gridcheck(*grid2d);
    checkIntersectionIterator(*grid2d);

    // dimworld == 3
    FoamGrid<3>* grid3d = make2Din3DHybridTestGrid<FoamGrid<3> >();

    gridcheck(*grid3d);
    checkIntersectionIterator(*grid3d);

    // dimworld == 3,  and a grid containing a T-Junction
    std::auto_ptr<FoamGrid<3> > gridTJunction( GmshReader<FoamGrid<3> >::read( dune_foamgrid_path + "tjunction-2d.msh", false, false ) );

    gridcheck(*gridTJunction);
    checkIntersectionIterator(*gridTJunction);
} 
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
