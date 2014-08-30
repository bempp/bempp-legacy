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

#include "hmat_parameters.hpp"

namespace Bempp {

ParameterEntry HMatParameters::defaultHMatAssemblyMode() {

  const std::string value("GlobalAssembly");
  return ParameterEntry(value, true, false,
                        "Specifies whether to assemble using global "
                        "or local dofs.");
}

ParameterEntry HMatParameters::defaultMinBlockSize() {

  const unsigned int value = 16;
  return ParameterEntry(value, true, false,
                        "Specifies the minimum block size below which "
                        "blocks are assumed to be dense.");
}

ParameterEntry HMatParameters::defaultMaxBlockSize() {

  const unsigned int value = 2048;
  return ParameterEntry(value, true, false,
                        "Specifies the maximum size of an admissible block.");
}

ParameterEntry HMatParameters::defaultEta() {

  constexpr double value = 1.2;
  return ParameterEntry(value, true, false,
                        "Specifies the block separation parameter eta.");
}

shared_ptr<ParameterList> HMatParameters::parameterList() {

  auto parameters = shared_ptr<ParameterList>(new ParameterList());

  parameters->setName("HMatParameters");
  parameters->setEntry("HMatAssemblyMode", defaultHMatAssemblyMode());
  parameters->setEntry("minBlockSize", defaultMinBlockSize());
  parameters->setEntry("maxBlockSize", defaultMaxBlockSize());
  parameters->setEntry("eta", defaultEta());

  return parameters;
}
}
