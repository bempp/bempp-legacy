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
        for pybasis, cybasis in dtypes.items()
        for pyresult, cyresult in dtypes.items()
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
% for pyname, ctype in dtypes.items():
        template<class T> std::string operator()(
            BoundaryOperator<${ctype}, T> const &_in) const {
            return "${pyname}";
        }
% endfor
    };
    struct ResultType: public boost::static_visitor<std::string> {
% for pyname, ctype in dtypes.items():
        template<class T> std::string operator()(
            BoundaryOperator<T, ${ctype}> const &_in) const {
            return "${pyname}";
        }
% endfor
    };

    struct Label: public boost::static_visitor<std::string> {
        template<class T> std::string operator()(T const &_in) const {
            return _in.label();
        }
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

% for op, opname in {'+': "Plus", '-': "Minus", '*': "Times"}.items():
    struct BoundaryOp${opname}BoundaryOp
        : boost::static_visitor<BoundaryOpVariants> {
          template<class BASIS, class RESULT> BoundaryOpVariants operator()(
                  BoundaryOperator<BASIS, RESULT> const &_a,
                  BoundaryOperator<BASIS, RESULT> const &_b) const {
              return BoundaryOpVariants(_a ${op} _b);
          }
          template<class BASIS0, class BASIS1, class RES0, class RES1>
              BoundaryOpVariants operator()(
                  BoundaryOperator<BASIS0, RES0> const &_a,
                  BoundaryOperator<BASIS1, RES1> const &_b) const {
                  throw std::runtime_error("Incompatible boundary operator");
              }
    };
% endfor
    template<class T> struct is_scalar : public boost::mpl::has_key<
        boost::mpl::set<
            float, double, std::complex<float>, std::complex<double> 
        >,
        T
    > {};
% for op, opname in {'/': "Divided", '*': "Times"}.items():
    template<class SCALAR> struct BoundaryOp${opname}Scalar
        : boost::static_visitor<BoundaryOpVariants> {
        SCALAR scalar;
        BoundaryOp${opname}Scalar(SCALAR _scalar) : scalar(_scalar) {}
        template<class BASIS, class RESULT>
            typename std::enable_if<
                std::is_convertible<SCALAR, BASIS>::value,
                BoundaryOpVariants
            > :: type
                operator()(BoundaryOperator<BASIS, RESULT> const &_a) const {
                    return BoundaryOpVariants(_a ${op} scalar);
                }
        template<class BASIS, class RESULT>
            typename std::enable_if<
                not std::is_convertible<SCALAR, BASIS>::value,
                BoundaryOpVariants
            > :: type
                operator()(BoundaryOperator<BASIS, RESULT> const &_a) const {
                    throw std::runtime_error("Incorrect scalar type");
                }
    };
% endfor


    public:
        BoundaryOpVariants() {}
        template<class BASIS, class RESULT>
            BoundaryOpVariants(BoundaryOperator<BASIS, RESULT> const &_bop)
                : operator_(_bop) {}

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
        std::string label() const {
            return boost::apply_visitor(Label(), operator_);
        }

% for op, opname in {'+': "Plus", '-': "Minus", '*': "Times"}.items():
        BoundaryOpVariants operator${op}(BoundaryOpVariants const &_a) const {
            return boost::apply_visitor(BoundaryOp${opname}BoundaryOp(),
                    operator_, _a.operator_);
        }
% endfor
% for op, opname in {'/': "Divided", '*': "Times"}.items():
    template<class T>
        typename std::enable_if<
            is_scalar<T> :: value,
            BoundaryOpVariants
        > :: type operator${op}(T const &_a) const {
            return boost::apply_visitor(BoundaryOp${opname}Scalar<T>(_a),
                    operator_);
        }
% endfor

    private:
        t_variant operator_;
};
}
#endif
