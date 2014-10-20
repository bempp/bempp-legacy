<% from space import dtypes, spaces %>
#ifndef BEMPP_PYTHON_SPACE_VARIANT_HPP
#define BEMPP_PYTHON_SPACE_VARIANT_HPP

#include "bempp/common/shared_ptr.hpp"
#include "types.h"
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

        bool isCompatible(SpaceVariants const &_other) const {
            return boost::apply_visitor(SpaceVariants::IsCompatible(),
                    space_, _other.space_);
        }
        bool isSame(SpaceVariants const &_other) const {
            return boost::apply_visitor(SpaceVariants::IsSame(),
                    space_, _other.space_);
        }

        t_variant & variants() { return space_; }
        t_variant const & variants() const { return space_; }
    private:
        t_variant space_;
};

#   undef BEMPP_EXPLICIT_CONSTRUCTOR

}
#endif
