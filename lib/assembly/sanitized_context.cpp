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

#include "sanitized_context.hpp"

#include "context.hpp"

#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <stdexcept>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Context<BasisFunctionType, ResultType>> sanitizedContext(
    shared_ptr<const Context<BasisFunctionType, ResultType>> context,
    bool localModeSupported, bool hybridModeSupported,
    const std::string &label) {
  if (!context)
    throw std::invalid_argument("sanitizedContext(): "
                                "context must not be null");

  const AssemblyOptions &assemblyOptions = context->assemblyOptions();
  const AcaOptions &acaOptions = assemblyOptions.acaOptions();
  if (assemblyOptions.assemblyMode() != AssemblyOptions::ACA)
    return context;
  bool switchToGlobalMode = false;
  if (acaOptions.mode == AcaOptions::LOCAL_ASSEMBLY && !localModeSupported) {
    if (acaOptions.reactionToUnsupportedMode == AcaOptions::ERROR)
      throw std::runtime_error("sanitizedContext(): "
                               "operator constructed by the function " +
                               label + " does not support assembly "
                                       "in the local ACA mode");
    else if (acaOptions.reactionToUnsupportedMode == AcaOptions::WARNING)
      std::cout << "sanitizedContext(): "
                   "operator constructed by the function " +
                       label +
                       " does not support assembly "
                       "in the local ACA mode; switching to the global ACA mode"
                << std::endl;
    switchToGlobalMode = true;
  } else if (acaOptions.mode == AcaOptions::HYBRID_ASSEMBLY &&
             !hybridModeSupported) {
    if (acaOptions.reactionToUnsupportedMode == AcaOptions::ERROR)
      throw std::runtime_error("sanitizedContext(): "
                               "operator constructed by the function " +
                               label + " does not support assembly "
                                       "in the hybrid ACA mode");
    else if (acaOptions.reactionToUnsupportedMode == AcaOptions::WARNING)
      std::cout
          << "sanitizedContext(): "
             "operator constructed by the function " +
                 label +
                 " does not support assembly "
                 "in the hybrid ACA mode; switching to the global ACA mode"
          << std::endl;
    switchToGlobalMode = true;
  }
  if (switchToGlobalMode) {
    // reset ACA mode to GLOBAL_ASSEMBLY
    AcaOptions newAcaOptions = acaOptions;
    newAcaOptions.mode = AcaOptions::GLOBAL_ASSEMBLY;
    AssemblyOptions newAssemblyOptions = assemblyOptions;
    newAssemblyOptions.switchToAcaMode(newAcaOptions);
    context.reset(new Context<BasisFunctionType, ResultType>(
        context->quadStrategy(), newAssemblyOptions));
  }
  return context;
}

#define INSTANTIATE_FUNCTION(BASIS, RESULT)                                    \
  template shared_ptr<const Context<BASIS, RESULT>> sanitizedContext(          \
      shared_ptr<const Context<BASIS, RESULT>>, bool, bool,                    \
      const std::string &)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FUNCTION);

} // namespace Bempp
