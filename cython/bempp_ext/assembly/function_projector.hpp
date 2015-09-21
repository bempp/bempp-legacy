#ifndef BEMPP_EXT_FUNCTION_PROJECTOR_HPP
#define BEMPP_EXT_FUNCTION_PROJECTOR_HPP

#include "bempp/fiber/surface_normal_and_domain_index_dependent_function.hpp"
#include "bempp/common/types.hpp"
#include <bempp/common/complex_aux.hpp>
#include "bempp/assembly/context.hpp"
#include "bempp/fiber/scalar_traits.hpp"
#include "bempp_ext/utils/py_types.hpp"
#include "bempp/assembly/local_assembler_construction_helper.hpp"
#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/basis.hpp"
#include "bempp/fiber/basis_data.hpp"
#include "bempp/fiber/collection_of_basis_transformations.hpp"
#include "bempp/fiber/function.hpp"
#include "bempp/fiber/local_assembler_for_grid_functions.hpp"
#include "bempp/fiber/quadrature_strategy.hpp"
#include "bempp/fiber/quadrature_strategy.hpp"
#include "bempp/fiber/raw_grid_geometry.hpp"
#include "bempp/grid/geometry_factory.hpp"
#include "bempp/grid/grid.hpp"
#include "bempp/grid/grid_view.hpp"
#include "bempp/grid/entity_iterator.hpp"
#include "bempp/grid/entity.hpp"
#include "bempp/grid/mapper.hpp"
#include "bempp/space/space.hpp"
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex>

namespace Bempp
{


template <typename ValueType_>
class PythonFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

    PythonFunctor(
        PyObject* callable,
        int argumentDimension, int resultDimension) :
            m_callable(callable),
            m_argumentDimension(argumentDimension),
            m_resultDimension(resultDimension)
            {

            Py_INCREF(m_callable);

            npy_intp pyArgumentDimension = m_argumentDimension;
            npy_intp pyResultDimension = m_resultDimension;

            m_x = PyArray_ZEROS(1,&pyArgumentDimension,NumpyType<CoordinateType>::value,1);
            m_normal = PyArray_ZEROS(1,&pyArgumentDimension,NumpyType<CoordinateType>::value,1);
            m_result = PyArray_ZEROS(1,&pyResultDimension,NumpyType<ValueType>::value,1);

            } 

    PythonFunctor(const PythonFunctor<ValueType>& other):
        m_argumentDimension(other.m_argumentDimension),
        m_resultDimension(other.m_resultDimension),m_x(other.m_x),m_normal(other.m_normal),
        m_result(other.m_result),m_callable(other.m_callable) {

            Py_INCREF(m_callable);
            Py_INCREF(m_x);
            Py_INCREF(m_normal);
            Py_INCREF(m_result);
            

        }

    ~PythonFunctor(){

        Py_DECREF(m_x);
        Py_DECREF(m_normal);
        Py_DECREF(m_result);
        Py_DECREF(m_callable);

    }


    int argumentDimension() const {
        return m_argumentDimension;
    }

    int resultDimension() const {
        return m_resultDimension;
    }

    void evaluate(const Eigen::Ref<Vector<CoordinateType>>& point, 
                  const Eigen::Ref<Vector<CoordinateType>>& normal,
                  int domainIndex, 
                  Eigen::Ref<Vector<ValueType>> result_) const
    {

        CoordinateType* xPtr = (CoordinateType*)PyArray_DATA(m_x);
        for (int i = 0; i< m_argumentDimension;++i) xPtr[i] = point(i);

        CoordinateType* normalPtr = (CoordinateType*)PyArray_DATA(m_normal);
        for (int i = 0; i< m_argumentDimension;++i) normalPtr[i] = normal(i);

        PyObject* argList = Py_BuildValue("OOiO",m_x,m_normal,domainIndex,m_result);
        PyObject* output = PyEval_CallObject(m_callable, argList);
        if (!output){
            Py_DECREF(argList);
            Py_DECREF(m_x);
            Py_DECREF(m_normal);
            Py_DECREF(m_result);
            Py_DECREF(m_callable);
            throw std::runtime_error("Error in evaluation of Python callable.");
        }
        Py_DECREF(output);
        Py_DECREF(argList);

        ValueType* resPtr = (ValueType*)PyArray_DATA(m_result);
        for (int i = 0; i< m_resultDimension;++i) result_(i) = resPtr[i];
    }

private:
    int m_argumentDimension;
    int m_resultDimension;
    PyObject* m_x;
    PyObject* m_normal;
    PyObject* m_result;
    PyObject* m_callable;

};


template <typename BasisFunctionType, typename ResultType>
Vector<ResultType> reallyCalculateProjections(
    const Space<BasisFunctionType> &dualSpace,
    Fiber::LocalAssemblerForGridFunctions<ResultType> &assembler,
    const AssemblyOptions &options) {
  // TODO: parallelise using TBB (the parameter options will then start be used)

  // Get the grid's leaf view so that we can iterate over elements
  const GridView &view = dualSpace.gridView();
  const size_t elementCount = view.entityCount(0);

  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> testGlobalDofs(elementCount);
  std::vector<std::vector<BasisFunctionType>> testLocalDofWeights(elementCount);

  // Gather global DOF lists
  const Mapper &mapper = view.elementMapper();
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    dualSpace.getGlobalDofs(element, testGlobalDofs[elementIndex],
                            testLocalDofWeights[elementIndex]);
    it->next();
  }

  // Make a vector of all element indices
  std::vector<int> testIndices(elementCount);
  for (size_t i = 0; i < elementCount; ++i)
    testIndices[i] = i;

  // Create the weak form's column vector
  Vector<ResultType> result(dualSpace.globalDofCount());
  result.setZero();

  std::vector<Vector<ResultType>> localResult;
  // Evaluate local weak forms
  assembler.evaluateLocalWeakForms(testIndices, localResult);

  // Loop over test indices
  for (size_t testIndex = 0; testIndex < elementCount; ++testIndex)
    // Add the integrals to appropriate entries in the global weak form
    for (size_t testDof = 0; testDof < testGlobalDofs[testIndex].size();
         ++testDof) {
      int testGlobalDof = testGlobalDofs[testIndex][testDof];
      if (testGlobalDof >= 0) // if it's negative, it means that this
                              // local dof is constrained (not used)
        result(testGlobalDof) +=
            conj(testLocalDofWeights[testIndex][testDof]) *
            localResult[testIndex](testDof);
    }

  // Return the vector of projections <phi_i, f>
  return result;
}

/** \brief Calculate projections of the function on the basis functions of
  the given dual space. */
template <typename BasisFunctionType, typename ResultType>
PyObject*
calculateProjections(const ParameterList& parameterList,
                     PyObject* callable,
                     const Space<BasisFunctionType> &dualSpace) {
  
  const Context<BasisFunctionType, ResultType> context(parameterList);


  auto globalFunction =  shared_ptr<Fiber::Function<ResultType>>(
     new Fiber::SurfaceNormalAndDomainIndexDependentFunction<PythonFunctor<ResultType>>(
         PythonFunctor<ResultType>(callable,3,dualSpace.codomainDimension())));

  const AssemblyOptions &options = context.assemblyOptions();

  // Prepare local assembler
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType> *>
      ShapesetPtrVector;
  typedef LocalAssemblerConstructionHelper Helper;

  shared_ptr<RawGridGeometry> rawGeometry;
  shared_ptr<GeometryFactory> geometryFactory;
  shared_ptr<Fiber::OpenClHandler> openClHandler;
  shared_ptr<ShapesetPtrVector> testShapesets;

  Helper::collectGridData(dualSpace, rawGeometry, geometryFactory);
  Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                            rawGeometry, openClHandler);
  Helper::collectShapesets(dualSpace, testShapesets);

  // Get reference to the test shapeset transformation
  const Fiber::CollectionOfShapesetTransformations<CoordinateType>
      &testTransformations = dualSpace.basisFunctionValue();

  typedef Fiber::LocalAssemblerForGridFunctions<ResultType> LocalAssembler;
  std::unique_ptr<LocalAssembler> assembler =
      context.quadStrategy()->makeAssemblerForGridFunctions(
          geometryFactory, rawGeometry, testShapesets,
          make_shared_from_ref(testTransformations),
          globalFunction, openClHandler);

  Vector<ResultType> result;
  result =  reallyCalculateProjections(dualSpace, *assembler, options);

  npy_intp pyResultDimension = dualSpace.globalDofCount();
  PyObject* pyResult = PyArray_ZEROS(1,&pyResultDimension,NumpyType<ResultType>::value,1);
  ResultType* resPtr = (ResultType*)PyArray_DATA(pyResult);

  for (int i = 0; i< pyResultDimension;++i) resPtr[i] = result(i);

  return pyResult;
}




} // namespace Bempp


#endif
