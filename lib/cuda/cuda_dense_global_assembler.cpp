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
#include "cuda_integrator.hpp"

#include "../assembly/discrete_dense_boundary_operator.hpp"

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/element_pair_topology.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/numerical_quadrature.hpp"

#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#include "../space/space.hpp"
#include "../grid/grid.hpp"

#include <tbb/tbb.h>
#include <algorithm>
#include <chrono>

namespace Bempp {

// Helper functions and classes
namespace {

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid() */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType> &space,
    std::vector<std::vector<GlobalDofIndex>> &globalDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) {

  // Get the grid's view so that we can iterate over elements
  const GridView &view = space.gridView();
  const int elementCount = view.entityCount(0);

  // Global DOF indices corresponding to local DOFs on elements
  globalDofs.clear();
  globalDofs.resize(elementCount);
  // Weights of the local DOFs on elements
  localDofWeights.clear();
  localDofWeights.resize(elementCount);

  // Gather global DOF lists
  const Mapper &mapper = view.elementMapper();
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    space.getGlobalDofs(element, globalDofs[elementIndex],
                        localDofWeights[elementIndex]);
    it->next();
  }
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
    std::vector<std::pair<
    std::pair<Matrix<typename ScalarTraits<BasisFunctionType>::RealType>,
    std::vector<typename ScalarTraits<BasisFunctionType>::RealType>>,
    std::pair<Matrix<typename ScalarTraits<BasisFunctionType>::RealType>,
    std::vector<typename ScalarTraits<BasisFunctionType>::RealType>>>>
    &regularQuadOrderCombinations,
    std::vector<std::vector<std::pair<
    const Fiber::Shapeset<BasisFunctionType>*,
    const Fiber::Shapeset<BasisFunctionType>*>>>
    &regularShapesetCombinations,
    std::vector<int> &singularElemPairTestIndices,
    std::vector<int> &singularElemPairTrialIndices) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef std::pair<const Shapeset*, const Shapeset*> ShapesetPair;
  typedef std::pair<Matrix<CoordinateType>, std::vector<CoordinateType>> QuadData;
  typedef std::pair<QuadData, QuadData> QuadDataPair;

  std::vector<const Shapeset*> testShapesets, trialShapesets;
  getAllShapesets(testSpace, testShapesets);
  if (&testSpace == &trialSpace) {
    trialShapesets = testShapesets;
  } else {
    getAllShapesets(trialSpace, trialShapesets);
  }

  const int testIndexCount = testIndices.size();
  const int trialIndexCount = trialIndices.size();

  unsigned int uniqueRegularQuadOrderCombinationCount = 0;
  std::vector<unsigned int> uniqueRegularShapesetCombinationCount;

  regularElemPairTestIndices.resize(1);
  regularElemPairTrialIndices.resize(1);
  regularElemPairTestIndices[0].resize(1);
  regularElemPairTrialIndices[0].resize(1);
  regularElemPairTestIndices[0][0].reserve(testIndexCount * trialIndexCount);
  regularElemPairTrialIndices[0][0].reserve(testIndexCount * trialIndexCount);

  regularQuadOrderCombinations.reserve(1);
  regularShapesetCombinations.resize(1);
  regularShapesetCombinations[0].reserve(1);

  singularElemPairTestIndices.reserve(testIndexCount * trialIndexCount);
  singularElemPairTrialIndices.reserve(testIndexCount * trialIndexCount);

  std::unique_ptr<GridView> testView = testSpace.grid()->leafView();
  std::unique_ptr<GridView> trialView = trialSpace.grid()->leafView();

  Fiber::RawGridGeometry<CoordinateType> testRawGeometry(
      testSpace.grid()->dim(), testSpace.grid()->dimWorld());
  Fiber::RawGridGeometry<CoordinateType> trialRawGeometry(
      trialSpace.grid()->dim(), trialSpace.grid()->dimWorld());

  testView->getRawElementData(
      testRawGeometry.vertices(), testRawGeometry.elementCornerIndices(),
      testRawGeometry.auxData(), testRawGeometry.domainIndices());
  trialView->getRawElementData(
      trialRawGeometry.vertices(), trialRawGeometry.elementCornerIndices(),
      trialRawGeometry.auxData(), trialRawGeometry.domainIndices());

  const Shapeset* testShapeset = testShapesets[testIndices[0]];
  const Shapeset* trialShapeset = trialShapesets[trialIndices[0]];
  ShapesetPair shapesetPair = std::make_pair(testShapeset, trialShapeset);
  const int testBasisOrder = testShapeset->order();
  const int trialBasisOrder = trialShapeset->order();
  const int testQuadOrder = testBasisOrder+4;
  const int trialQuadOrder = trialBasisOrder+4;
  QuadDataPair quadDataPair;
  Fiber::fillSingleQuadraturePointsAndWeights(
      3, testQuadOrder, quadDataPair.first.first, quadDataPair.first.second);
  Fiber::fillSingleQuadraturePointsAndWeights(
      3, trialQuadOrder, quadDataPair.second.first, quadDataPair.second.second);
  regularQuadOrderCombinations.push_back(quadDataPair);
  regularShapesetCombinations[0].push_back(shapesetPair);

  // Sort element pairs regarding regular and singular integrals
  for (int testIndex = 0; testIndex < testIndexCount; ++testIndex) {
    for (int trialIndex = 0; trialIndex < trialIndexCount; ++trialIndex) {

      // Get corner indices of the specified elements
      Vector<int> testElementCornerIndices =
          testRawGeometry.elementCornerIndices(testIndices[testIndex]);
      Vector<int> trialElementCornerIndices =
          trialRawGeometry.elementCornerIndices(trialIndices[trialIndex]);

      Fiber::ElementPairTopology topology;
      if (testSpace.grid().get() == trialSpace.grid().get()) {

        topology = Fiber::determineElementPairTopologyIn3D(
            testElementCornerIndices,
            trialElementCornerIndices);

      } else {

        topology.testVertexCount = testElementCornerIndices.rows();
        topology.trialVertexCount = trialElementCornerIndices.rows();
        topology.type = Fiber::ElementPairTopology::Disjoint;
      }

      if (topology.type == Fiber::ElementPairTopology::Disjoint) {

        // Element pair related to regular integrals
        regularElemPairTestIndices[0][0].push_back(testIndices[testIndex]);
        regularElemPairTrialIndices[0][0].push_back(trialIndices[trialIndex]);

      } else {

        // Element pair related to singular integrals
        singularElemPairTestIndices.push_back(testIndices[testIndex]);
        singularElemPairTrialIndices.push_back(trialIndices[trialIndex]);
      }
    }
  }
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &assembler,
    const Context<BasisFunctionType, ResultType> &context) {

  std::cout << "Hello, this is CudaDenseGlobalAssembler::assembleDetachedWeakForm()!" << std::endl;

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
  trialIndices.reserve(trialElementCount);

  for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex) {
    const int trialDofCount = trialGlobalDofs[trialIndex].size();
    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
      if (trialGlobalDof >= 0) {
        trialIndices.push_back(trialIndex);
        break;
      }
    }
  }

//  std::cout << "testIndices = " << std::endl;
//  for (int i = 0; i < testIndices.size(); ++i) {
//    std::cout << testIndices[i] << " " << std::flush;
//  }
//  std::cout << std::endl;
//
//  std::cout << "trialIndices = " << std::endl;
//  for (int i = 0; i < trialIndices.size(); ++i) {
//    std::cout << trialIndices[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  // Create the operator's matrix
  Matrix<ResultType> result(testSpace.globalDofCount(),
                            trialSpace.globalDofCount());
  result.setZero();

  // Get element pair index vectors and regular integration data
  std::vector<std::vector<std::vector<int>>> regularElemPairTestIndices;
  std::vector<std::vector<std::vector<int>>> regularElemPairTrialIndices;
  std::vector<QuadDataPair> regularQuadOrderCombinations;
  std::vector<std::vector<ShapesetPair>> regularShapesetCombinations;
  std::vector<int> singularElemPairTestIndices;
  std::vector<int> singularElemPairTrialIndices;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  getSortedElementPairs(
      testSpace, trialSpace, testIndices, trialIndices,
      regularElemPairTestIndices, regularElemPairTrialIndices,
      regularQuadOrderCombinations, regularShapesetCombinations,
      singularElemPairTestIndices, singularElemPairTrialIndices);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for getSortedElementPairs() = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

//  std::cout << "number of singular element pairs: " << singularElemPairTestIndices.size() << std::endl;
//  std::cout << "number of regular element pairs: " << regularElemPairTestIndices[0][0].size() << std::endl;

//  std::cout << "singularElemPairIndices = " << std::endl;
//  for (int i = 0; i < singularElemPairTestIndices.size(); ++i) {
//    std::cout << singularElemPairTestIndices[i] << " "
//              << singularElemPairTrialIndices[i] << std::endl;
//  }
//  std::cout << std::endl;
//
//  std::cout << "regularElemPairIndices = " << std::endl;
//  for (int i = 0; i < regularElemPairTestIndices[0][0].size(); ++i) {
//    std::cout << regularElemPairTestIndices[0][0][i] << " "
//              << regularElemPairTrialIndices[0][0][i] << std::endl;
//  }
//  std::cout << std::endl;

//  std::cout << "testQuadPoints = " << std::endl;
//  std::cout << regularQuadOrderCombinations[0].first.first << std::endl;
//  std::cout << "testQuadWeights = " << std::endl;
//  for (int i = 0; i < regularQuadOrderCombinations[0].first.second.size(); ++i) {
//    std::cout << regularQuadOrderCombinations[0].first.second[i] << " " << std::flush;
//  }
//  std::cout << std::endl;
//  std::cout << "trialQuadPoints = " << std::endl;
//  std::cout << regularQuadOrderCombinations[0].second.first << std::endl;
//  std::cout << "trialQuadWeights = " << std::endl;
//  for (int i = 0; i < regularQuadOrderCombinations[0].second.second.size(); ++i) {
//    std::cout << regularQuadOrderCombinations[0].second.second[i] << " " << std::flush;
//  }
//  std::cout << std::endl;
//
//  std::cout << "number of test shape functions: " << regularShapesetCombinations[0][0].first->size() << std::endl;
//  std::cout << "number of trial shape functions: " << regularShapesetCombinations[0][0].second->size() << std::endl;

  begin = std::chrono::steady_clock::now();

  // Evaluate singular integrals over selected element pairs
  std::vector<Matrix<ResultType>> singularResult;
  assembler.evaluateLocalWeakForms(singularElemPairTestIndices,
                                   singularElemPairTrialIndices,
                                   singularResult);

  end = std::chrono::steady_clock::now();
  std::cout << "Time for evaluateSingularWeakForms() = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  std::vector<std::vector<std::vector<Matrix<ResultType>>>> regularResult;

  // Push raw grid data to the device
  shared_ptr<CudaGrid> testGrid = testSpace.grid()->pushToDevice(0);
  shared_ptr<CudaGrid> trialGrid = trialSpace.grid()->pushToDevice(0);

  const unsigned int quadOrderCombinationCount = regularQuadOrderCombinations.size();

  regularResult.resize(quadOrderCombinationCount);

  begin = std::chrono::steady_clock::now();

  for (int quadOrderCombination = 0; quadOrderCombination < quadOrderCombinationCount; ++quadOrderCombination) {

    const unsigned int shapesetCombinationCount = regularShapesetCombinations[quadOrderCombination].size();

    Matrix<CoordinateType> localTestQuadPoints =
        regularQuadOrderCombinations[quadOrderCombination].first.first;
    Matrix<CoordinateType> localTrialQuadPoints =
        regularQuadOrderCombinations[quadOrderCombination].second.first;
    std::vector<CoordinateType> testQuadWeights =
        regularQuadOrderCombinations[quadOrderCombination].first.second;
    std::vector<CoordinateType> trialQuadWeights =
        regularQuadOrderCombinations[quadOrderCombination].second.second;

    regularResult[quadOrderCombination].resize(shapesetCombinationCount);

    for (int shapesetCombination = 0; shapesetCombination < shapesetCombinationCount; ++shapesetCombination) {

      const unsigned int elemPairCount = regularElemPairTestIndices[quadOrderCombination][shapesetCombination].size();

      const Shapeset &testShapeset = *(regularShapesetCombinations[quadOrderCombination][shapesetCombination].first);
      const Shapeset &trialShapeset = *(regularShapesetCombinations[quadOrderCombination][shapesetCombination].second);

      Fiber::CudaIntegrator<BasisFunctionType, ResultType> cudaIntegrator(
          localTestQuadPoints,
          localTrialQuadPoints,
          testQuadWeights,
          trialQuadWeights,
          testShapeset,
          trialShapeset,
          testGrid,
          trialGrid);

      regularResult[quadOrderCombination][shapesetCombination].resize(elemPairCount);

      // Create chunks of element pairs according to a maximum number of
      // element pairs active on the device
      const unsigned int maxActiveElemPairCount = 100e06;

      std::vector<int> elemPairChunkTestIndices;
      std::vector<int> elemPairChunkTrialIndices;
      elemPairChunkTestIndices.reserve(maxActiveElemPairCount);
      elemPairChunkTrialIndices.reserve(maxActiveElemPairCount);

      std::vector<Matrix<ResultType>*> regularResultChunk;
      regularResultChunk.reserve(maxActiveElemPairCount);

      unsigned int activeElemPairCount = 0;

      for (int elemPair = 0; elemPair < elemPairCount; ++elemPair) {

//        std::cout << "elemPair = " << elemPair << std::endl;

        const int elemPairTestIndex =
            regularElemPairTestIndices[quadOrderCombination][shapesetCombination][elemPair];
        const int elemPairTrialIndex =
            regularElemPairTrialIndices[quadOrderCombination][shapesetCombination][elemPair];

        if (activeElemPairCount < maxActiveElemPairCount) {

          elemPairChunkTestIndices.push_back(elemPairTestIndex);
          elemPairChunkTrialIndices.push_back(elemPairTrialIndex);
          regularResultChunk.push_back(&regularResult[quadOrderCombination][shapesetCombination][elemPair]);
          activeElemPairCount++;
        }

//        std::cout << "elemPairChunkIndices = " << std::endl;
//        for (int i = 0; i < elemPairChunkTestIndices.size(); ++i) {
//          std::cout << elemPairChunkTestIndices[i] << " "
//                    << elemPairChunkTrialIndices[i] << std::endl;
//        }

        if (activeElemPairCount == maxActiveElemPairCount || elemPair == elemPairCount-1) {

          end = std::chrono::steady_clock::now();
          std::cout << "Time for creating chunk = "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                    << " ms" << std::endl;

          begin = std::chrono::steady_clock::now();

          // Evaluate regular integrals over selected element pairs
          cudaIntegrator.integrate(elemPairChunkTestIndices,
                                   elemPairChunkTrialIndices,
                                   regularResultChunk);

          end = std::chrono::steady_clock::now();
          std::cout << "Time for integrate() = "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                    << " ms" << std::endl;

//          std::cout << "regularResultChunk = " << std::endl;
//          for (int i = 0; i < elemPairChunkTestIndices.size(); ++i) {
//            std::cout << *(regularResultChunk[i]) << std::endl;
//          }

          elemPairChunkTestIndices.clear();
          elemPairChunkTrialIndices.clear();
          regularResultChunk.clear();
          activeElemPairCount = 0;

          begin = std::chrono::steady_clock::now();
        }
      }
    }
  }

  // Global assembly
  {
    begin = std::chrono::steady_clock::now();

    // Loop over singular element pairs
    for (int singularElemPair = 0; singularElemPair < singularElemPairTestIndices.size(); ++singularElemPair) {

      const int testIndex = singularElemPairTestIndices[singularElemPair];
      const int trialIndex = singularElemPairTrialIndices[singularElemPair];

      const int testDofCount = testGlobalDofs[testIndex].size();
      const int trialDofCount = trialGlobalDofs[trialIndex].size();

      // Add the integrals to appropriate entries in the operator's matrix
      for (int testDof = 0; testDof < testDofCount; ++testDof) {

        int testGlobalDof = testGlobalDofs[testIndex][testDof];
        if (testGlobalDof < 0)
          continue;

        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {

          int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
          if (trialGlobalDof < 0)
            continue;

          assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);
          assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);

          result(testGlobalDof, trialGlobalDof) +=
            conj(testLocalDofWeights[testIndex][testDof]) *
            trialLocalDofWeights[trialIndex][trialDof] *
            singularResult[singularElemPair](testDof, trialDof);
        }
      }
    }

    // Loop over regular element pairs
    for (int quadOrderCombination = 0; quadOrderCombination < quadOrderCombinationCount; ++quadOrderCombination) {

      const unsigned int shapesetCombinationCount = regularResult[quadOrderCombination].size();
      for (int shapesetCombination = 0; shapesetCombination < shapesetCombinationCount; ++shapesetCombination) {

        const unsigned int elemPairCount = regularResult[quadOrderCombination][shapesetCombination].size();
        for (int elemPair = 0; elemPair < elemPairCount; ++elemPair) {

          const int testIndex = regularElemPairTestIndices[quadOrderCombination][shapesetCombination][elemPair];
          const int trialIndex = regularElemPairTrialIndices[quadOrderCombination][shapesetCombination][elemPair];

          const int testDofCount = testGlobalDofs[testIndex].size();
          const int trialDofCount = trialGlobalDofs[trialIndex].size();

          // Add the integrals to appropriate entries in the operator's matrix
          for (int testDof = 0; testDof < testDofCount; ++testDof) {

            int testGlobalDof = testGlobalDofs[testIndex][testDof];
            if (testGlobalDof < 0)
              continue;

            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {

              int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
              if (trialGlobalDof < 0)
                continue;

              assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);
              assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);

              result(testGlobalDof, trialGlobalDof) +=
                conj(testLocalDofWeights[testIndex][testDof]) *
                trialLocalDofWeights[trialIndex][trialDof] *
                regularResult[quadOrderCombination][shapesetCombination][elemPair](testDof, trialDof);
            }
          }
        }
      }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Time for global assembly = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << " ms" << std::endl;
  }

  std::cout << "result (cudadense) = " << std::endl;
  std::cout << result << std::endl;

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
