#include "standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/standard_local_assembler_for_grid_functions_on_surfaces_imp.hpp"
#include "../fiber/standard_local_assembler_for_identity_operator_on_surface_imp.hpp"
#include "../fiber/standard_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

namespace Bempp
{

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
        StandardLocalAssemblerFactoryForOperatorsOnSurfaces);

} // namespace Bempp
