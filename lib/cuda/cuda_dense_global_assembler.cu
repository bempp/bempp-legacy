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

#include "../assembly/discrete_dense_boundary_operator.hpp"
#include "../common/types.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

//// TODO Helper functions
//namespace {
//
///** Build a list of lists of global DOF indices corresponding to the local DOFs
// *  on each element of space.grid() */
//template <typename BasisFunctionType>
//void gatherGlobalDofs(
//    const Space<BasisFunctionType> &space,
//    std::vector<std::vector<GlobalDofIndex>> &globalDofs,
//    std::vector<std::vector<BasisFunctionType>> &localDofWeights) {
//
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
//}
//
//} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &assembler,
    const Context<BasisFunctionType, ResultType> &context) {

  // Push the test and trial grid data to the device
  shared_ptr<CudaGrid> testGrid = testSpace.grid().pushToDevice();
  shared_ptr<CudaGrid> trialGrid = trialSpace.grid().pushToDevice();
  testGrid->setupGeometry();
  trialGrid->setupGeometry();

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

  // Create the operator's matrix
  Matrix<ResultType> result(testSpace.globalDofCount(),
                            trialSpace.globalDofCount());
  result.setZero();

  std::vector<Matrix<ResultType>> localResult;
  for (int trialIndex = 0; trialIndex != trialElementCount; ++trialIndex) {
    // Handle this trial element only if it contributes to any global DOFs
    bool skipTrialElement = true;
    const int trialDofCount = trialGlobalDofs[trialIndex].size();
    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
      if (trialGlobalDof >= 0) {
        skipTrialElement = false;
        break;
      }
    }
    if (skipTrialElement)
      continue;

    // Evaluate integrals over pairs of the current trial element and
    // all the test elements
    assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices, trialIndex,
                                     ALL_DOFS, localResult);

    // Global assembly
    {
      // Loop over test indices
      for (int row = 0; row < testElementCount; ++row) {
        const int testIndex = testIndices[row];
        const int testDofCount = testGlobalDofs[testIndex].size();
        // Add the integrals to appropriate entries in the operator's matrix
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
          int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
          if (trialGlobalDof < 0)
            continue;
          for (int testDof = 0; testDof < testDofCount; ++testDof) {
            int testGlobalDof = testGlobalDofs[testIndex][testDof];
            if (testGlobalDof < 0)
              continue;
            assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);
            assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);
            result(testGlobalDof, trialGlobalDof) +=
              conj(testLocalDofWeights[testIndex][testDof]) *
              trialLocalDofWeights[trialIndex][trialDof] *
              localResult[row](testDof, trialDof);
          }
        }
      }
    }
  }

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
