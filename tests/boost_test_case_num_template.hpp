// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

//  This file contains fragments of the source code of the Boost::Test library
// (specifically, test_case_template.hpp and boost/test/unit_test_suite.hpp).

//  This is the original licence statement accompanying these files:

//  (C) Copyright Gennadiy Rozental 2003-2005.
//  Distributed under the Boost Software License, Version 1.0.
//  (http://www.boost.org/LICENSE_1_0.txt)

/** \file

    Adds support for test case templates instantiated with a
    sequence of test types, such as boost::mpl::int_<...>, whose member
    "value" defines an integer constant. This value is then used
    instead of the type id to provide a human-readable test name.

    \example
    #include "boost_test_case_num_template.hpp"
    #include <boost/mpl/int.hpp>
    typedef boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<1>, boost::mpl::int_<2> >
    my_type_list;

    BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(my_test, T, my_type_list) {
        // inside the test, the value of the parameter (0, 1 or 2) is available from T::value.
        BOOST_CHECK(T::value == 0 || T::value == 1 || T::value == 2);
        BOOST_CHECK_EQUAL(T::value, 3); // will fail
    }

    A typical message output by Boost::Test if such a test fails will be
    filename.cpp(266): error in "my_test<0>": check T::value == 3 failed [0 != 3]
*/
// ***************************************************************************

#ifndef BOOST_TEST_TEST_CASE_NUM_TEMPLATE_HPP
#define BOOST_TEST_TEST_CASE_NUM_TEMPLATE_HPP

// Boost.Test
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_case_template.hpp>

// Boost
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/type.hpp>
#include <boost/type_traits/is_const.hpp>

// STL
#include <sstream>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

#define BOOST_TEST_CASE_NUM_TEMPLATE( name, typelist )                          \
    boost::unit_test::ut_detail::num_template_test_case_gen<name,typelist >(    \
        BOOST_TEST_STRINGIZE( name ) )                                      \
/**/

//____________________________________________________________________________//

// A numeric variant of BOOST_TEST_CASE_TEMPLATE_FUNCTION is not
// currently provided.

//____________________________________________________________________________//

namespace boost {

namespace unit_test {

namespace ut_detail {

//____________________________________________________________________________//

// ************************************************************************** //
// ************           generate_test_case_4_num_type          ************ //
// ************************************************************************** //

template<typename Generator,typename TestCaseTemplate>
struct generate_test_case_4_num_type {
    explicit    generate_test_case_4_num_type( const_string tc_name, Generator& G )
    : m_test_case_name( tc_name )
    , m_holder( G )
    {}

    template<typename TestType>
    void        operator()( mpl::identity<TestType> )
    {
        std::string full_name;
        assign_op( full_name, m_test_case_name, 0 );

        // This is the essential fragment changed wrt. boost_test_case_template.hpp
        std::stringstream sstream;
        sstream << full_name;
        sstream << '<';
        sstream << TestType::value;
        sstream << '>';
        full_name = sstream.str();

        m_holder.m_test_cases.push_back( 
            new test_case( full_name, test_case_template_invoker<TestCaseTemplate,TestType>() ) );
    }

private:
    // Data members
    const_string    m_test_case_name;
    Generator&      m_holder;
};

// ************************************************************************** //
// ************              test_case_num_template              ************ //
// ************************************************************************** //

template<typename TestCaseTemplate,typename TestTypesList>
class num_template_test_case_gen : public test_unit_generator {
public:
    // Constructor
    num_template_test_case_gen( const_string tc_name )
    {
        typedef generate_test_case_4_num_type<num_template_test_case_gen<TestCaseTemplate,TestTypesList>,
                                          TestCaseTemplate
        > single_test_gen;
        mpl::for_each<TestTypesList,mpl::make_identity<mpl::_> >( single_test_gen( tc_name, *this ) );
    }

    test_unit* next() const
    {
        if( m_test_cases.empty() )
            return 0;
    
        test_unit* res = m_test_cases.front();
        m_test_cases.pop_front();

        return res;
    }

    // Data members
    mutable std::list<test_unit*> m_test_cases;
};

//____________________________________________________________________________//

} // namespace ut_detail

} // unit_test

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

//____________________________________________________________________________//

// ************************************************************************** //
// ************        BOOST_AUTO_TEST_CASE_NUM_TEMPLATE         ************ //
// ************************************************************************** //

#define BOOST_AUTO_TEST_CASE_NUM_TEMPLATE( test_name, type_name, TL )   \
template<typename type_name>                                            \
struct test_name : public BOOST_AUTO_TEST_CASE_FIXTURE                  \
{ void test_method(); };                                                \
                                                                        \
struct BOOST_AUTO_TC_INVOKER( test_name ) {                             \
    template<typename TestType>                                         \
    static void run( boost::type<TestType>* = 0 )                       \
    {                                                                   \
        test_name<TestType> t;                                          \
        t.test_method();                                                \
    }                                                                   \
};                                                                      \
                                                                        \
BOOST_AUTO_TU_REGISTRAR( test_name )(                                   \
    boost::unit_test::ut_detail::num_template_test_case_gen<            \
        BOOST_AUTO_TC_INVOKER( test_name ),TL >(                        \
          BOOST_STRINGIZE( test_name ) ) );                             \
                                                                        \
template<typename type_name>                                            \
void test_name<type_name>::test_method()                                \
/**/

#endif // BOOST_TEST_TEST_CASE_NUM_TEMPLATE_HPP

