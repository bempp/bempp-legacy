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

#ifndef bempp_assembly_options_hpp
#define bempp_assembly_options_hpp

#include "../common/common.hpp"

#include <string>

#include "../fiber/opencl_options.hpp"
#include "../fiber/parallelization_options.hpp"

namespace Bempp
{

/** \ingroup assembly
 *  \brief Adaptive cross approximation (ACA) parameters.
 */
struct AcaOptions
{
    /** \brief Initialize ACA parameters to default values. */
    AcaOptions();

    /** \brief Estimate of the desired approximation accuracy.
     *
     *  Default value: 1e-4. */
    double eps;
    /** \brief Cluster-pair admissibility parameter.
     *
     *  Default value: 1.2. */
    double eta;
    /** \brief Minimum size of blocks approximated with ACA.
     *
     *  Matrix blocks whose smaller side is less than this value will be stored
     *  in the dense format (ACA will not be attempted on these blocks).
     *
     *  Default value: 16. */
    unsigned int minimumBlockSize;
    /** \brief Maximum allowed block size.
     *
     *  Matrix blocks with side larger than this value will be split. This can
     *  be used to limit the time devoted to (serial) processing of an
     *  individual block, especially when a large number of threads is used.
     *
     *  Default value: UINT_MAX. */
    unsigned int maximumBlockSize;
    /** \brief Maximum rank of blocks stored in the low-rank format.
     *
     *  Blocks judged to have higher rank will be stored in the dense format.
     *
     *  Default value: UINT_MAX. */
    unsigned int maximumRank;
    /** \brief Do global assembly before ACA?
     *
     *  If true, ACA is performed on the assembled operator (i.e. the matrix is
     *  indexed with global degrees of freedom). Otherwise ACA is performed on
     *  matrix indexed with local degrees of freedom, and the global assembly
     *  is done by pre- and postmultiplying the H matrix with sparse matrices
     *  mapping local DOFs to global DOFs and vice versa.
     *
     *  Default value: true. */
    bool globalAssemblyBeforeCompression;
    /** \brief Recompress ACA matrix after construction?
     *
     *  If true, blocks of H matrices are agglomerated in an attempt to reduce
     *  memory consumption.
     *
     *  Warning: this procedure is not parallelised yet, therefore it may be
     *  slow.
     *
     *  Default value: false. */
    bool recompress;
    /** \brief If true, hierarchical matrix structure will be written in
     *  PostScript format at the end of the assembly procedure.
     *
     *  Default value: false. */
    bool outputPostscript;
    /** \brief Name of the output PostScript file.
     *
     *  \see outputPostscript.
     *
     *  Default value: "aca.ps". */
    std::string outputFname;
    /** \brief Estimate of the magnitude of typical entries of the matrix to be
     *  approximated.
     *
     *  Usually does not need to be changed. Default value: 1. */
    double scaling;
};

using Fiber::OpenClOptions;
using Fiber::ParallelizationOptions;

/** \ingroup assembly
 *  \brief Options determining how weak-form assembly is done.
 */
class AssemblyOptions
{
public:
    AssemblyOptions();

    /** @name Operator representation
      @{ */

    enum Mode {
        AUTO = -1,
        NO = 0,
        YES = 1
    };

    enum Representation {
        DENSE, ACA, FMM, SPARSE
    };

    void switchToDense();
    void switchToAca(const AcaOptions& acaOptions);
    void switchToFmm();

    Representation operatorRepresentation() const {
        return m_representation;
    }

    const AcaOptions& acaOptions() const {
        return m_acaOptions;
    }

    /** @}
      @name Parallelization
      @{ */

    void switchToOpenCl(const OpenClOptions& openClOptions);
    void switchToTbb(int maxThreadCount = AUTO);

    const ParallelizationOptions& parallelizationOptions() const {
        return m_parallelizationOptions;
    }

    /** @}
      @name Others
      @{ */

    /** \brief Are singular integrals cached?

      Possible settings:
        * YES: precalculate singular integrals
        * NO: do not precalculate singular integrals
        * AUTO: implementation-defined
          (currently singular integrals are cached in ACA mode only)
      */
    void setSingularIntegralCaching(Mode mode);

    Mode singularIntegralCaching() const {
        return m_singularIntegralCaching;
    }

    /** @} */

private:
    Representation m_representation;
    AcaOptions m_acaOptions;
    ParallelizationOptions m_parallelizationOptions;
    Mode m_singularIntegralCaching;
};

} // namespace Bempp

#endif
