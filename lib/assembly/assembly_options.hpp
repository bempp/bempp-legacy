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

#include "../fiber/opencl_options.hpp"

namespace Bempp
{

struct AcaOptions
{
    double eps;
    double eta;
    int minimumBlockSize;
    int maximumRank;
    bool recompress;
};

using Fiber::OpenClOptions;

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
    void switchToSparse();

    Representation operatorRepresentation() const {
        return m_representation;
    }

    const AcaOptions& acaOptions() const {
        return m_acaOptions;
    }

    /** @}
      @name Parallelism
      @{ */

    enum Parallelism {
        TBB_AND_OPEN_MP, OPEN_CL
    };

    void switchToOpenCl(const OpenClOptions& openClOptions);
    void switchToTbbAndOpenMp(int maxThreadCount = AUTO);

    Parallelism parallelism() const {
        return m_parallelism;
    }

    const OpenClOptions& openClOptions() const {
        return m_openClOptions;
    }

    int maxThreadCount() const {
        return m_maxThreadCount;
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
    Parallelism m_parallelism;
    OpenClOptions m_openClOptions;
    int m_maxThreadCount;
    Mode m_singularIntegralCaching;
};

} // namespace Bempp

#endif
