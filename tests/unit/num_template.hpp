// Copyright (C) 2011-2012 by the BEM++ Authors
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

#ifndef bempp_num_template_hpp
#define bempp_num_template_hpp

#include "boost_test_case_num_template.hpp"
#include <boost/mpl/list.hpp>
#include <boost/mpl/int.hpp>

typedef boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<1>, boost::mpl::int_<2> >
list_0_to_2;

typedef boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<1>, boost::mpl::int_<2>, boost::mpl::int_<3> >
list_0_to_3;

typedef boost::mpl::list<boost::mpl::int_<1>, boost::mpl::int_<2> > list_1_to_2;

typedef boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<1> > list_0_to_1;

typedef boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<3> > list_0_and_3;

#endif
