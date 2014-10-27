<%
from bempp_operators import bops, compatible_dtypes, dtypes
%>\
#ifndef BEMPP_PYTHON_ASSEMBLY_PYTHON_H
#define BEMPP_PYTHON_ASSEMBLY_PYTHON_H

#include <boost/variant/multivisitors.hpp>

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/assembly/numerical_quadrature_strategy.hpp"
#include "bempp/assembly/context.hpp"
#include "bempp/assembly/variants.hpp"
#include "bempp/space/types.h"

% for description in bops.values():
#include "${description['header']}"
% endfor

namespace {
    using namespace Bempp;

    //! Basis and Result are compatible
    template<class BASIS, class RESULT> struct is_compatible
        : std::integral_constant<bool, false> {};
% for basis, results in compatible_dtypes.items():
%   for result in results:
    template<> struct is_compatible<${dtypes[basis]}, ${dtypes[result]}>
        : std::integral_constant<bool, true> {};
%   endfor
% endfor


% for opname, description in bops.items():
%     if description['implementation'] == 'standard':
    // \brief Factory functor for boundary operator
    // \details Works with a BoundaryOpVariants (and hence boost::variants) to 
    // create a boundary operator via the factory function ${opname}
    template<class RESULT> struct CreateBoundaryOps${opname}
        : public boost::static_visitor<BoundaryOpVariants> {
            CreateBoundaryOps${opname}(
                const Fiber::AccuracyOptionsEx& _accuracyOptions,
                const AssemblyOptions& _assemblyOptions,
                const std::string& _label,
                int _symmetry
            ) : boost::static_visitor<BoundaryOpVariants>(), // needed by intel
                accuracyOptions(_accuracyOptions),
                assemblyOptions(_assemblyOptions),
                label(_label), symmetry(_symmetry) {}

        template<class SPACE> typename std::enable_if<
            // Make sure instantiate only if basis and result types are
            // compatible
            is_compatible<typename SPACE::BasisFunctionType, RESULT> :: value,
            BoundaryOpVariants
        > :: type operator()(
            const shared_ptr<SPACE>& _domain,
            const shared_ptr<SPACE>& _range,
            const shared_ptr<SPACE>& _dualToRange
          ) const {
            typedef typename SPACE::BasisFunctionType t_Basis;
            typedef Context<t_Basis, RESULT> t_Context;
            typedef typename t_Context::QuadratureStrategy t_Strategy;
            typedef NumericalQuadratureStrategy<t_Basis, RESULT> t_NumStrategy;

            const shared_ptr<const t_Strategy> quadStrategy(
                new const t_NumStrategy(accuracyOptions)
            );
            const shared_ptr<const t_Context> context(
                new const t_Context(quadStrategy, assemblyOptions)
            );
            // At least clang can't figure constructing shared_ptr to const
            // Space automatically. So do it by hand.
            shared_ptr<Space<t_Basis> const> domain(_domain);
            shared_ptr<Space<t_Basis> const> range(_range);
            shared_ptr<Space<t_Basis> const> dualToRange(_dualToRange);
            BoundaryOpVariants result;
            result.set(Bempp::${opname}(
                context,
                domain, range, dualToRange,
                label, symmetry
            ));
            return result;
          }

        template<class A, class B, class C> BoundaryOpVariants operator()(
            const shared_ptr<A>& _domain,
            const shared_ptr<B>& _range,
            const shared_ptr<C>& _dualToRange
        ) const {
            throw std::runtime_error("domain, range and/or dual to range "
                    "spaces are incompatible");
        }
        template<class BASIS> typename std::enable_if<
            not is_compatible<BASIS, RESULT> :: value,
            BoundaryOpVariants
        > :: type operator()(
            const shared_ptr<Space<BASIS> const>& _domain,
            const shared_ptr<Space<BASIS> const>& _range,
            const shared_ptr<Space<BASIS> const>& _dualToRange
        ) const {
            throw std::runtime_error("Incompatible basis and result types");
        }

        Fiber::AccuracyOptionsEx accuracyOptions;
        AssemblyOptions assemblyOptions;
        std::string label;
        int symmetry;
    };
    //! Creates a ${opname} operator and stores is in a BoundaryOpVariants
    template<class RESULT>
      BoundaryOpVariants ${description['c_creator']}(
            const Fiber::AccuracyOptionsEx& accuracyOptions,
            const AssemblyOptions& assemblyOptions,
            const SpaceVariants& domain,
            const SpaceVariants& range,
            const SpaceVariants& dualToRange,
            const std::string& label,
            int symmetry
        ) {
            return boost::apply_visitor(
                CreateBoundaryOps${opname}<RESULT>(accuracyOptions,
                    assemblyOptions, label, symmetry),
                domain.variants(), range.variants(), dualToRange.variants()
            );
        }
%     endif
% endfor
}

#endif
