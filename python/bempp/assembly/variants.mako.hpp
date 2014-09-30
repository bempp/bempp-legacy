<% from bempp_operators import dtypes, compatible_dtypes, bops %>
#ifndef BEMPP_PYTHON_ASSEMBLY_VARIANTS_HPP
#define BEMPP_PYTHON_ASSEMBLY_VARIANTS_HPP

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/space/types.h"
#include "bempp/space/variant.hpp"
#include <boost/variant.hpp>
#include <type_traits>

namespace Bempp {

//! A shared pointer that can hold const and non-const versions
class BoundaryOpVariants {
<%
    variants = [(cybasis, cyresult)
        for pybasis, cybasis in dtypes.iteritems()
        for pyresult, cyresult in dtypes.iteritems()
        if pyresult in compatible_dtypes[pybasis]
    ]
%>
    typedef boost::variant<
% for cybasis, cyresult in variants:
        BoundaryOperator<${cybasis}, ${cyresult}>
            ${',' if not loop.last else ''}
% endfor
    > t_variant;

    struct BasisType: public boost::static_visitor<std::string> {
% for pyname, ctype in dtypes.iteritems():
        template<class T> std::string operator()(
            BoundaryOperator<${ctype}, T> const &_in) const {
            return "${pyname}";
        }
% endfor
    };
    struct ResultType: public boost::static_visitor<std::string> {
% for pyname, ctype in dtypes.iteritems():
        template<class T> std::string operator()(
            BoundaryOperator<T, ${ctype}> const &_in) const {
            return "${pyname}";
        }
% endfor
    };


% for space in ['Domain', 'Range', 'DualToRange']:
    struct ${space}: public boost::static_visitor<SpaceVariants> {
        template<class BASIS, class RESULT> SpaceVariants operator()(
                BoundaryOperator<BASIS, RESULT> const &_in) const {
            <% name = space[0].lower() + space[1:] %>
            shared_ptr<Space<BASIS> const> result = _in.${name}();
            if(not result)
                throw std::runtime_error(
                        "Operator does not have a space defined");
            SpaceVariants space;
            space.set(result);
            return space;
        }
    };
%endfor

    public:
        BoundaryOpVariants() {}

        template<class Basis, class Result>
            void set(BoundaryOperator<Basis, Result> const &_in) { operator_ = _in; }
        template<class Basis, class Result>
            void set() { operator_ = BoundaryOperator<Basis, Result>(); }

        std::string basisType() const {
            return boost::apply_visitor(BasisType(), operator_);
        }
        std::string resultType() const {
            return boost::apply_visitor(ResultType(), operator_);
        }
        SpaceVariants domain() const {
            return boost::apply_visitor(Domain(), operator_);
        }
        SpaceVariants range() const {
            return boost::apply_visitor(Range(), operator_);
        }
        SpaceVariants dualToRange() const {
            return boost::apply_visitor(DualToRange(), operator_);
        }

    private:
        t_variant operator_;
};

}
#endif
