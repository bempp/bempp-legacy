<% 
from data_types import dtypes
from space import spaces 

%>
#ifndef BEMPP_PYTHON_SPACE_VARIANTS_HPP
#define BEMPP_PYTHON_SPACE_VARIANTS_HPP

#include "bempp/common/shared_ptr.hpp"
#include "bempp/space/space.hpp"
#include "bempp/utils/py_types.hpp"
#include <boost/variant.hpp>
#include <type_traits>
#include <type_traits>

namespace Bempp {

# ifdef __INTEL_COMPILER
#   define BEMPP_EXPLICIT_CONSTRUCTOR(NAME, TYPE) \
        NAME() : boost::static_visitor<TYPE>() {}
# else
#   define BEMPP_EXPLICIT_CONSTRUCTOR(NAME, TYPE) 
# endif

//! A shared pointer that can hold const and non-const versions
class SpaceVariants {
    typedef boost::variant<
% for ctype in dtypes.values():
        shared_ptr< Space<${ctype}> >,
        shared_ptr< const Space<${ctype}> > ${',' if not loop.last else ''}
% endfor
    > t_variant;

    struct GetGrid : public boost::static_visitor< shared_ptr<Grid const> > {
        BEMPP_EXPLICIT_CONSTRUCTOR(GetGrid, shared_ptr<Grid const>);
        template<class T> shared_ptr<Grid const>
            operator()( shared_ptr<T> const &_in) const { return _in->grid(); }
    };

    struct GetCodomainDimension : public boost::static_visitor< int > {
        BEMPP_EXPLICIT_CONSTRUCTOR(GetCodomainDimension, int);
        template<typename T> int operator()(shared_ptr<T> const &_in) const { return _in->codomainDimension();}
    };


    struct IsDiscontinuous : public boost::static_visitor< bool > {
        BEMPP_EXPLICIT_CONSTRUCTOR(IsDiscontinuous, bool);
        template<typename T> bool operator()(shared_ptr<T> const &_in) const { return _in->isDiscontinuous();}
    };


    struct GetDomainDimension : public boost::static_visitor< int > {
        BEMPP_EXPLICIT_CONSTRUCTOR(GetDomainDimension, int);
        template<typename T> int operator()(shared_ptr<T> const &_in) const { return _in->domainDimension();}
    };


    struct GetGlobalDofCount : public boost::static_visitor< size_t > {
        BEMPP_EXPLICIT_CONSTRUCTOR(GetGlobalDofCount, size_t);
        template<typename T> size_t operator()(shared_ptr<T> const &_in) const { return _in->globalDofCount();}
    };

    struct GetFlatLocalDofCount : public boost::static_visitor< size_t > {
        BEMPP_EXPLICIT_CONSTRUCTOR(GetFlatLocalDofCount, size_t);
        template<typename T> size_t operator()(shared_ptr<T> const &_in) const { return _in->flatLocalDofCount();}
    };


    struct IsCompatible : public boost::static_visitor<bool> {
        BEMPP_EXPLICIT_CONSTRUCTOR(IsCompatible, bool);
        template<class T0, class T1>
            typename std::enable_if<
                std::is_same<
                    typename std::remove_const<T0>::type,
                    typename std::remove_const<T1>::type
                > :: value, bool
            > :: type operator()(shared_ptr<T0> const &_this,
                shared_ptr<T1> const &_that) const {
            return _this->spaceIsCompatible(*_that);
        }
        template<class T0, class T1>
            typename std::enable_if<
                not std::is_same<
                    typename std::remove_const<T0>::type,
                    typename std::remove_const<T1>::type
                > :: value, bool
            > :: type operator()(shared_ptr<T0> const &_this,
                shared_ptr<T1> const &_that) const {
            return false;
        }
    };

    struct IsSame : public boost::static_visitor<bool> {
        BEMPP_EXPLICIT_CONSTRUCTOR(IsSame, bool);
        template<class T0, class T1>
            typename std::enable_if<
                std::is_same<
                    typename std::remove_const<T0>::type,
                    typename std::remove_const<T1>::type
                > :: value, bool
            > :: type operator()(shared_ptr<T0> const &_this,
                shared_ptr<T1> const &_that) const {
            return _this.get() == _that.get();
        }
        template<class T0, class T1>
            typename std::enable_if<
                not std::is_same<
                    typename std::remove_const<T0>::type,
                    typename std::remove_const<T1>::type
                > :: value, bool
            > :: type operator()(shared_ptr<T0> const &_this,
                shared_ptr<T1> const &_that) const {
            return false;
        }
    };

    struct DType: public boost::static_visitor<std::string> {
        BEMPP_EXPLICIT_CONSTRUCTOR(DType, std::string);
% for pyname, ctype in dtypes.items():
        std::string operator()(
                shared_ptr<Space<${ctype}> const> const &_in) const {
            return "${pyname}";
        }
% endfor
    };

    template <typename BasisFunctionType>
    friend shared_ptr<const Space<BasisFunctionType>> _py_get_space_ptr(const SpaceVariants& space_variant);

    public:
        SpaceVariants() {}

        template<class T>
            void set(shared_ptr<T> const &_in) { space_ = _in; }
        template<class T>
            void set(shared_ptr<T const> const &_in) { space_ = _in; }

        std::string dtype() const {
            return boost::apply_visitor(SpaceVariants::DType(), space_);
        }

        shared_ptr<Grid const> grid() const {
            return boost::apply_visitor(SpaceVariants::GetGrid(), space_);
        }

        int codomainDimension() const {
            return boost::apply_visitor(SpaceVariants::GetCodomainDimension(), space_);
        }

        bool isCompatible(SpaceVariants const &_other) const {
            return boost::apply_visitor(SpaceVariants::IsCompatible(),
                    space_, _other.space_);
        }
        bool isSame(SpaceVariants const &_other) const {
            return boost::apply_visitor(SpaceVariants::IsSame(),
                    space_, _other.space_);
        }

        int domainDimension() const {
            return boost::apply_visitor(SpaceVariants::GetDomainDimension(),space_);
        }

        size_t globalDofCount() const {
            return boost::apply_visitor(SpaceVariants::GetGlobalDofCount(),space_);
        }

        size_t flatLocalDofCount() const {
            return boost::apply_visitor(SpaceVariants::GetFlatLocalDofCount(),space_);
        }

        t_variant & variants() { return space_; }
        t_variant const & variants() const { return space_; }
    private:
        t_variant space_;
};

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>> _py_get_space_ptr(const SpaceVariants& space_variant)
{
    shared_ptr<const Space<BasisFunctionType>> res;
    bool success = true;

    try {
        res = boost::get<shared_ptr<Space<BasisFunctionType>>>(space_variant.space_);
    }
    catch (...) {
        success = false;
    }
    if (!success)
        res = boost::get<shared_ptr<const Space<BasisFunctionType>>>(space_variant.space_);
    return res;
}

template<typename BasisFunctionType>
arma::Mat<typename Fiber::ScalarTraits<BasisFunctionType>::RealType> _py_space_get_global_dof_interp_points(
        const SpaceVariants& space_variant){

    shared_ptr<const Space<BasisFunctionType>> space_ptr = _py_get_space_ptr<BasisFunctionType>(space_variant);
    arma::Mat<typename Fiber::ScalarTraits<BasisFunctionType>::RealType> result;
    space_ptr->getGlobalDofInterpolationPoints(result);
    return result;
}


template<typename BasisFunctionType>
arma::Mat<typename Fiber::ScalarTraits<BasisFunctionType>::RealType> _py_space_get_global_dof_normals(
        const SpaceVariants& space_variant){

    shared_ptr<const Space<BasisFunctionType>> space_ptr = _py_get_space_ptr<BasisFunctionType>(space_variant);
    arma::Mat<typename Fiber::ScalarTraits<BasisFunctionType>::RealType> result;
    space_ptr->getNormalsAtGlobalDofInterpolationPoints(result);
    return result;
}

#   undef BEMPP_EXPLICIT_CONSTRUCTOR

}
#endif
