#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <dune/common/exceptions.hh>
#include <Python.h>
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
      catch (...)
      {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
      }
    }
}
