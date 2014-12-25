#ifndef BEMPP_PYTHON_UTILS_HPP
#define BEMPP_PYTHON_UTILS_HPP

#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <dune/common/exceptions.hh>
#include <Python.h>
#include <Teuchos_ParameterList.hpp>
#include <algorithm>
#include <sstream>
#include <complex.h>

namespace Bempp {
    inline static void catch_exception() {
      try {
        // latest Python exn passes through, ignore the current one
        if (not PyErr_Occurred()) throw;
      } catch (const Dune::IOError& exn) {
        PyErr_SetString(PyExc_IOError, exn.what().c_str());
      } catch (const Dune::Exception& exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what().c_str());
      } catch (const std::bad_alloc& exn) {
        PyErr_SetString(PyExc_MemoryError, exn.what());
      } catch (const std::bad_cast& exn) {
        PyErr_SetString(PyExc_TypeError, exn.what());
      } catch (const std::domain_error& exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
      } catch (const std::invalid_argument& exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
      } catch (const std::ios_base::failure& exn) {
        PyErr_SetString(PyExc_IOError, exn.what());
      } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
      } catch (const std::overflow_error& exn) {
        PyErr_SetString(PyExc_OverflowError, exn.what());
      } catch (const std::range_error& exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
      } catch (const std::underflow_error& exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
      } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what());
      }
      catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
      }
    }

    inline static int parameter_list_length(const Teuchos::ParameterList& parameters)
    {
        int i = 0;
        for (auto it = std::begin(parameters); it!= std::end(parameters);++it) ++i;
        return i;
    }

    inline static std::vector<std::string> parameter_names(const Teuchos::ParameterList& parameters)
    {
        std::vector<std::string> names;
        for (auto it = std::begin(parameters);it!=std::end(parameters);++it)
        {
            names.push_back(it->first);
        }
        return names;
    }

    inline static std::string print_parameters(const Teuchos::ParameterList& parameters)
    {

        std::stringstream ss;
        ss << parameters;
        return ss.str();
    }

}
#endif

