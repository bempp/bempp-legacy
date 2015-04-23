<% from bempp_operators import dtypes, compatible_dtypes, bops %>
#ifndef BEMPP_PYTHON_ASSEMBLY_VARIANTS_HPP
#define BEMPP_PYTHON_ASSEMBLY_VARIANTS_HPP

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/utils/py_types.hpp"
#include "bempp/space/py_space_variants.hpp"
#include "bempp/assembly/discrete_sparse_boundary_operator.hpp"
#include "bempp/utils/py_types.hpp"
#include "fiber/scalar_traits.hpp"
#include <Python.h>
#include "numpy/arrayobject.h"
#include <boost/variant.hpp>
#include <type_traits>

namespace Bempp {

# ifdef __INTEL_COMPILER
#   define BEMPP_EXPLICIT_CONSTRUCTOR(NAME, TYPE) \
        NAME() : boost::static_visitor<TYPE>() {}
# else
#   define BEMPP_EXPLICIT_CONSTRUCTOR(NAME, TYPE) 
# endif


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
        BEMPP_EXPLICIT_CONSTRUCTOR(BasisType, std::string);
% for pyname, ctype in dtypes.items():
        template<class T> std::string operator()(
            BoundaryOperator<${ctype}, T> const &_in) const {
            return "${pyname}";
        }
% endfor
    };
    struct ResultType: public boost::static_visitor<std::string> {
        BEMPP_EXPLICIT_CONSTRUCTOR(ResultType, std::string);
% for pyname, ctype in dtypes.items():
        template<class T> std::string operator()(
            BoundaryOperator<T, ${ctype}> const &_in) const {
            return "${pyname}";
        }
% endfor
    };

    struct Label: public boost::static_visitor<std::string> {
        BEMPP_EXPLICIT_CONSTRUCTOR(Label, std::string);
        template<class T> std::string operator()(T const &_in) const {
            return _in.label();
        }
    };


% for space in ['Domain', 'Range', 'DualToRange']:
    struct Valid${space}: public boost::static_visitor<bool> {
        BEMPP_EXPLICIT_CONSTRUCTOR(Valid${space}, bool);
        template<class BASIS, class RESULT> bool operator()(
                BoundaryOperator<BASIS, RESULT> const &_in) const {
            <% name = space[0].lower() + space[1:] %>
            return (bool) _in.${name}();
        }
    };

    struct ${space}: public boost::static_visitor<SpaceVariants> {
        BEMPP_EXPLICIT_CONSTRUCTOR(${space}, SpaceVariants);
        template<class BASIS, class RESULT> SpaceVariants operator()(
                BoundaryOperator<BASIS, RESULT> const &_in) const {
            <% name = space[0].lower() + space[1:] %>
            shared_ptr<Space<BASIS> const> result = _in.${name}();
            SpaceVariants space;
            if(result) space.set(result);
            return space;
        }
    };
%endfor

% for op, opname in {'+': "Plus", '-': "Minus", '*': "Times"}.items():
    struct BoundaryOp${opname}BoundaryOp
        : boost::static_visitor<BoundaryOpVariants> {
          BEMPP_EXPLICIT_CONSTRUCTOR(
              BoundaryOp${opname}BoundaryOp, BoundaryOpVariants);
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
        BoundaryOp${opname}Scalar(SCALAR _scalar)
            : boost::static_visitor<BoundaryOpVariants>(), scalar(_scalar) {}
        template<class BASIS, class RESULT>
            typename std::enable_if<
                std::is_convertible<SCALAR, RESULT>::value,
                BoundaryOpVariants
            > :: type
                operator()(BoundaryOperator<BASIS, RESULT> const &_a) const {
                    return BoundaryOpVariants(_a ${op} scalar);
                }
        template<class BASIS, class RESULT>
            typename std::enable_if<
                not std::is_convertible<SCALAR, RESULT>::value,
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
% for space in ['Domain', 'Range', 'DualToRange']:
        <% name = space[0].lower() + space[1:] %>
	SpaceVariants ${name}() const {
            return boost::apply_visitor(${space}(), operator_);
        }
	bool valid${space}() const {
            return boost::apply_visitor(Valid${space}(), operator_);
        }
% endfor
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

        template<typename BasisFunctionType,typename ResultType>
        friend shared_ptr<const DiscreteBoundaryOperator<ResultType>>
        boundary_op_variant_weak_form(const BoundaryOpVariants& variant);

        t_variant operator_;
};
#   undef BEMPP_EXPLICIT_CONSTRUCTOR

template<typename BasisFunctionType,typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>> 
boundary_op_variant_weak_form(const BoundaryOpVariants& variant)
{

    return boost::get<BoundaryOperator<BasisFunctionType,ResultType>>(
            variant.operator_).weakForm();
}


}
#endif
