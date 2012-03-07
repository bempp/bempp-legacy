#ifndef bempp_stl_io_hpp
#define bempp_stl_io_hpp

#include <vector>
#include <iostream>

template <typename T>
std::ostream& operator<< (std::ostream& dest, const std::vector<T>& v)
{
    for (int i = 0; i < v.size(); ++i)
        dest << v[i] << '\n';
    dest.flush();
}

#endif
