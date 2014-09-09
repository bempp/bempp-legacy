<%
from bempp_operators import bops
%>\
#ifndef BEMPP_PYTHON_ASSEMBLY_PYTHON_H
#define BEMPP_PYTHON_ASSEMBLY_PYTHON_H

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/assembly/numerical_quadrature_strategy.hpp"
#include "bempp/assembly/context.hpp"

% for description in bops.itervalues():
#include "${description['header']}"
% endfor

namespace {
    using namespace Bempp;

    // Constructs C++ object in python-managed memory
    template<class BASIS, class RESULT>
        inline void inplace_boundary_operator(void *_memory) {
            new(_memory) BoundaryOperator<BASIS, RESULT>();
        }

    // Deconstructs C++ object
    template<class BASIS, class RESULT>
        inline void deconstructor(void *_memory) {
            static_cast<BoundaryOperator<BASIS, RESULT>*>(_memory)
                ->~BoundaryOperator<BASIS, RESULT>();
        }

% for opname, description in bops.iteritems():
%     if description['implementation'] == 'standard':
    template<class BASIS, class RESULT> void ${description['c_creator']}(
            BoundaryOperator<BASIS, RESULT>& _output,
            const Fiber::AccuracyOptionsEx& accuracyOptions,
            const AssemblyOptions& assemblyOptions,
            const shared_ptr<const Space<BASIS> >& domain,
            const shared_ptr<const Space<BASIS> >& range,
            const shared_ptr<const Space<BASIS> >& dualToRange,
            const std::string& label,
            int symmetry
        ) {
            typedef Context<BASIS, RESULT> t_Context;
            typedef typename t_Context::QuadratureStrategy t_Strategy;
            typedef NumericalQuadratureStrategy<BASIS, RESULT> t_NumStrategy;

            const shared_ptr<const t_Strategy> quadStrategy(
                new const t_NumStrategy(accuracyOptions)
            );
            const shared_ptr<const t_Context> context(
                new const t_Context(quadStrategy, assemblyOptions)
            );
            _output = Bempp::${opname}(
                context,
                domain, range, dualToRange,
                label, symmetry
            );
        }
%     endif
% endfor
}
#endif
