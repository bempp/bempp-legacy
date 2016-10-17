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

#include "cuda_dense_global_assembler.hpp"

#include "cuda_grid.hpp"
//#include "cuda_integrator.hpp"

#include "../assembly/discrete_dense_boundary_operator.hpp"

#include "../common/types.hpp"
//#include "../common/eigen_support.hpp"
//
#include "../fiber/explicit_instantiation.hpp"
//#include "../fiber/element_pair_topology.hpp"
//#include "../fiber/raw_grid_geometry.hpp"
//#include "../fiber/element_pair_topology.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
//#include "../fiber/shapeset.hpp"
//#include "../fiber/numerical_quadrature.hpp"
//
//#include "../grid/grid.hpp" // CAUSES COMPILER ERROR
//#include "../grid/grid_view.hpp" // CAUSES COMPILER ERROR

#include "../space/space.hpp"

namespace Bempp {

// Helper functions and classes
namespace {

/** TODO Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid(). */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType> &space,
    std::vector<std::vector<GlobalDofIndex>> &globalDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) {

//  // Get the grid's view so that we can iterate over elements
//  const GridView &view = space.gridView();
//  const int elementCount = view.entityCount(0);
//
//  // Global DOF indices corresponding to local DOFs on elements
//  globalDofs.clear();
//  globalDofs.resize(elementCount);
//  // Weights of the local DOFs on elements
//  localDofWeights.clear();
//  localDofWeights.resize(elementCount);
//
//  // Gather global DOF lists
//  const Mapper &mapper = view.elementMapper();
//  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
//  while (!it->finished()) {
//    const Entity<0> &element = it->entity();
//    const int elementIndex = mapper.entityIndex(element);
//    space.getGlobalDofs(element, globalDofs[elementIndex],
//                        localDofWeights[elementIndex]);
//    it->next();
//  }
}

// Get sorted element pair index vectors regarding regular and singular integrals
// as well as quadrature order and shapeset combinations for regular ones
template <typename BasisFunctionType>
void getSortedElementPairs(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    std::vector<int> &testIndices,
    std::vector<int> &trialIndices,
    std::vector<std::vector<std::vector<int>>> &regularElemPairTestIndices,
    std::vector<std::vector<std::vector<int>>> &regularElemPairTrialIndices,
    std::vector<std::pair<std::pair<Matrix<double>, std::vector<double>>, std::pair<Matrix<double>, std::vector<double>>>>
    &regularQuadOrderCombinations,
    std::vector<std::vector<std::pair<Fiber::Shapeset<BasisFunctionType>*, Fiber::Shapeset<BasisFunctionType>*>>>
    &regularShapesetCombinations,
    std::vector<int> &singularElemPairTestIndices,
    std::vector<int> &singularElemPairTrialIndices) {

  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef std::pair<Shapeset*, Shapeset*> ShapesetPair;
  typedef std::pair<Matrix<double>, std::vector<double>> QuadData;
  typedef std::pair<QuadData, QuadData> QuadDataPair;

  std::vector<const Shapeset*> testShapesets, trialShapesets;
  getAllShapesets(testSpace, testShapesets);
  if (&testSpace == &trialSpace) {
    trialShapesets = testShapesets;
  } else {
    getAllShapesets(trialSpace, trialShapesets);
  }

//  const int testIndexCount = testIndices.size();
//  const int trialIndexCount = trialIndices.size();
//
//  unsigned int uniqueRegularShapesetCombinationCount = 0;
//
//  regularElemPairTestIndices.reserve(testIndexCount * trialIndexCount);
//  regularElemPairTrialIndices.reserve(testIndexCount * trialIndexCount);
//  singularElemPairTestIndices.reserve(testIndexCount * trialIndexCount);
//  singularElemPairTrialIndices.reserve(testIndexCount * trialIndexCount);
//
//  std::unique_ptr<GridView> testView = testSpace.grid()->leafView();
//  std::unique_ptr<GridView> trialView = trialSpace.grid()->leafView();
//
//  Fiber::RawGridGeometry<double> testRawGeometry(
//      testSpace.grid()->dim(), testSpace.grid()->dimWorld());
//  Fiber::RawGridGeometry<double> trialRawGeometry(
//      trialSpace.grid()->dim(), trialSpace.grid()->dimWorld());
//
//  testView->getRawElementData(
//      testRawGeometry.vertices(), testRawGeometry.elementCornerIndices(),
//      testRawGeometry.auxData(), testRawGeometry.domainIndices());
//  trialView->getRawElementData(
//      trialRawGeometry.vertices(), trialRawGeometry.elementCornerIndices(),
//      trialRawGeometry.auxData(), trialRawGeometry.domainIndices());
//
//  // Sort element pairs regarding regular and singular integrals
//  for (int testIndex = 0; testIndex < testIndexCount; ++testIndex) {
//    for (int trialIndex = 0; trialIndex < trialIndexCount; ++trialIndex) {
//
//      // Get corner indices of the specified elements
//      Vector<int> testElementCornerIndices =
//          testRawGeometry.elementCornerIndices(testIndices[testIndex]);
//      Vector<int> trialElementCornerIndices =
//          trialRawGeometry.elementCornerIndices(trialIndices[trialIndex]);
//
//      Fiber::ElementPairTopology topology;
//      if (testSpace.grid().get() == trialSpace.grid().get()) {
//
//        topology = Fiber::determineElementPairTopologyIn3D(
//            testElementCornerIndices,
//            trialElementCornerIndices);
//
//      } else {
//
//        topology.testVertexCount = testElementCornerIndices.rows();
//        topology.trialVertexCount = trialElementCornerIndices.rows();
//        topology.type = Fiber::ElementPairTopology::Disjoint;
//      }
//
//      if (topology.type == Fiber::ElementPairTopology::Disjoint) {
//
//        // Element pair related to regular integrals
//        Shapeset* testShapeset = testShapesets[testIndex];
//        Shapeset* trialShapeset = trialShapesets[trialIndex];
//
//        int testBasisOrder = testShapeset->order();
//        int trialBasisOrder = trialShapeset->order();
//        testQuadOrder = testBasisOrder;
//        trialQuadOrder = trialBasisOrder;
//
//        bool putDown = false;
//        for (int i = 0; i < uniqueRegularShapesetCombinationCount; ++i) {
//
//          if (testShapeset == (*m_testShapesets)[regularElemPairTestIndices[i][0]] &&
//              trialShapeset == (*m_trialShapesets)[regularElemPairTrialIndices[i][0]]) {
//
//            regularElemPairTestIndices[i].push_back(testIndices[testIndex]);
//            regularElemPairTrialIndices[i].push_back(trialIndices[trialIndex]);
//
//            putDown = true;
//          }
//        }
//        if (putDown == false) {
//
//          regularElemPairTestIndices.resize(++uniqueRegularShapesetCombinationCount);
//          regularElemPairTrialIndices.resize(++uniqueRegularShapesetCombinationCount);
//
//          regularElemPairTestIndices[uniqueRegularShapesetCombinationCount].push_back(testIndex);
//          regularElemPairTrialIndices[uniqueRegularShapesetCombinationCount].push_back(trialIndex);
//
//          fillSingleQuadraturePointsAndWeights(
//              topology.testVertexCount, testQuadOrder, testQuadPoints, testQuadWeights);
//          fillSingleQuadraturePointsAndWeights(
//              topology.trialVertexCount, trialQuadOrder, trialQuadPoints, trialQuadWeights);
//        }
//
//      } else {
//
//        // Element pair related to singular integrals
//        singularElemPairTestIndices.push_back(testIndices[testIndex]);
//        singularElemPairTrialIndices.push_back(trialIndices[trialIndex]);
//      }
//    }
//  }
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &assembler,
    const Context<BasisFunctionType, ResultType> &context) {

  std::cout << "assemble detached weak form in CUDA dense mode" << std::endl;

  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> testGlobalDofs, trialGlobalDofs;
  std::vector<std::vector<BasisFunctionType>> testLocalDofWeights,
    trialLocalDofWeights;

  gatherGlobalDofs(testSpace, testGlobalDofs, testLocalDofWeights);

  if (&testSpace == &trialSpace) {
    trialGlobalDofs = testGlobalDofs;
    trialLocalDofWeights = testLocalDofWeights;
  } else {
    gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);
  }

  const int testElementCount = testGlobalDofs.size();
  const int trialElementCount = trialGlobalDofs.size();

  // Enumerate the test elements that contribute to at least one global DOF
  std::vector<int> testIndices;
  testIndices.reserve(testElementCount);

  for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
    const int testDofCount = testGlobalDofs[testIndex].size();
    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      int testGlobalDof = testGlobalDofs[testIndex][testDof];
      if (testGlobalDof >= 0) {
        testIndices.push_back(testIndex);
        break;
      }
    }
  }

  // Enumerate the trial elements that contribute to at least one global DOF
  std::vector<int> trialIndices;
  testIndices.reserve(trialElementCount);

  for (int trialIndex = 0; trialIndex < testElementCount; ++trialIndex) {
    const int trialDofCount = trialGlobalDofs[trialIndex].size();
    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
      if (trialGlobalDof >= 0) {
        trialIndices.push_back(trialIndex);
        break;
      }
    }
  }

  // Create the operator's matrix
  Matrix<ResultType> result(testSpace.globalDofCount(),
                            trialSpace.globalDofCount());
  result.setZero();

  // Get element pair index vectors
  std::vector<std::vector<std::vector<int>>> regularElemPairTestIndices;
  std::vector<std::vector<std::vector<int>>> regularElemPairTrialIndices;
  std::vector<QuadDataPair> regularQuadOrderCombinations;
  std::vector<std::vector<ShapesetPair>> regularShapesetCombinations;

  std::vector<int> singularElemPairTestIndices;
  std::vector<int> singularElemPairTrialIndices;

  getSortedElementPairs(
      testSpace, trialSpace, testIndices, trialIndices,
      regularElemPairTestIndices, regularElemPairTrialIndices,
      regularQuadOrderCombinations, regularShapesetCombinations,
      singularElemPairTestIndices, singularElemPairTrialIndices);

  // Evaluate singular integrals over selected element pairs
  std::vector<Matrix<ResultType>> singularResult;
  assembler.evaluateLocalWeakForms(singularElemPairTestIndices,
                                   singularElemPairTrialIndices,
                                   singularResult);

//  // Push raw grid data to the device
//  shared_ptr<CudaGrid> testGrid = testSpace.grid()->pushToDevice(0);
//  shared_ptr<CudaGrid> trialGrid = trialSpace.grid()->pushToDevice(0);

  const unsigned int quadOrderCombinationCount = regularElemPairTestIndices.size();

  for (int quadOrderCombination = 0; quadOrderCombination < quadOrderCombinationCount; ++quadOrderCombination) {

    const unsigned int shapesetCombinationCount = regularElemPairTestIndices[quadOrderCombination].size();

    Matrix<double> localTestQuadPoints;
    Matrix<double> localTrialQuadPoints;
    std::vector<double> testQuadWeights;
    std::vector<double> trialQuadWeights;

//    Fiber::CudaIntegrator<BasisFunctionType, ResultType> cudaIntegrator(
//        localTestQuadPoints,
//        localTrialQuadPoints,
//        testQuadWeights,
//        trialQuadWeights,
//        testGrid,
//        trialGrid);

    for (int shapesetCombination = 0; shapesetCombination < shapesetCombinationCount; ++shapesetCombination) {

      const unsigned int elemPairCount = regularElemPairTestIndices[quadOrderCombination][shapesetCombination].size();

      // Create chunks of element pairs according to a maximum number of element
      // data on the device
      const unsigned int maxDeviceElemCount = 1000;

      for (int elemPair = 0; elemPair < elemPairCount; ++elemPair) {

        unsigned int deviceElemCount = 0;

        std::vector<int> testDeviceElemIndices;
        std::vector<int> trialDeviceElemIndices;
        testDeviceElemIndices.reserve(maxDeviceElemCount);
        trialDeviceElemIndices.reserve(maxDeviceElemCount);

        std::vector<int> elemPairChunkTestIndices;
        std::vector<int> elemPairChunkTrialIndices;
        elemPairChunkTestIndices.reserve(elemPairCount);
        elemPairChunkTrialIndices.reserve(elemPairCount);

        const int elemPairTestIndex =
            regularElemPairTestIndices[quadOrderCombination][shapesetCombination][elemPair];
        const int elemPairTrialIndex =
            regularElemPairTrialIndices[quadOrderCombination][shapesetCombination][elemPair];

        // TODO
        if (true) {

          elemPairChunkTestIndices.push_back(elemPairTestIndex);
          elemPairChunkTrialIndices.push_back(elemPairTrialIndex);
        }

        // TODO
        if (elemPairChunkTestIndices.size() == elemPairCount) {

//          // Setup geometry data for selected elements on the device
//          testGrid->setupElements(testDeviceElemIndices);
//          trialGrid->setupElements(trialDeviceElemIndices);
//
//          // Evaluate regular integrals over selected element pairs
//          cudaIntegrator.integrate(elemPairChunkTestIndices,
//                                   elemPairChunkTrialIndices,
//                                   testShapeset, trialShapeset,
//                                   regularResultChunk);
//
//          // Free element data on the device
//          testGrid->freeElementData();
//          trialGrid->freeElementData();

          // TODO Assemble regular results
        }
      }
    }
  }

    // TODO Global assembly
//    {
//      // Loop over test indices
//      for (int row = 0; row < testElementCount; ++row) {
//        const int testIndex = testIndices[row];
//        const int testDofCount = testGlobalDofs[testIndex].size();
//        // Add the integrals to appropriate entries in the operator's matrix
//        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
//          int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
//          if (trialGlobalDof < 0)
//            continue;
//          for (int testDof = 0; testDof < testDofCount; ++testDof) {
//            int testGlobalDof = testGlobalDofs[testIndex][testDof];
//            if (testGlobalDof < 0)
//              continue;
//            assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);
//            assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);
//            result(testGlobalDof, trialGlobalDof) +=
//              conj(testLocalDofWeights[testIndex][testDof]) *
//              trialLocalDofWeights[trialIndex][trialDof] *
//              localResult[row](testDof, trialDof);
//          }
//        }
//      }
//    }
//  }

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
