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

#include "global_parameters.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace Bempp {

ParameterList GlobalParameters::parameterList() {

  ParameterList parameters;

  parameters.set(
      "maxThreadCount", static_cast<int>(-1),
      "(int) Specifies the maximum number of threads to use "
      "-1 denotes automatic selection of number of threads via TBB.");

  parameters.set("boundaryOperatorAssemblyType", std::string("dense"),
                  "(string) Default assembly type for boundary operators. "
                  "Allowed values are dense and hmat.");

  parameters.set("potentialOperatorAssemblyType", std::string("dense"),
          "(string) Default assembly type for potential oeprators. "
          "Allowed values are dense and hmat.");

  parameters.set("verbosityLevel",
          static_cast<int>(0),
          "(int) Default Verbosity of BEM++. Supported values are "
          "-5 (low verbosity), 0 (default), 5 (high verbosity)");

  parameters.set("enableSingularIntegralCaching", 
          true,
          "(bool) If true then singular integrals are pre-calculated and cached "
          "before the boundary operator assembly");


  parameters.set("enableBlasInQuadrature",
          std::string("auto"),
          "(std::string) Specifies whether to use BLAS in quadrature routines. "
          " If set to auto use only for quadratic and higher order basis functions "
          "(default). If set to  no disable and if set to yes  always enable.");

   
  ParameterList& quadratureOrders = parameters.sublist("QuadratureOrders");

  quadratureOrders.set("quadratureOrdersAreRelative",
         static_cast<bool>(true),
        "(bool) Specify whether quadrature options are relative to default "
        "orders. ");

  quadratureOrders.set("doubleSingular",static_cast<int>(0),
          "(int) Order for singular double integrals.");

  auto createQuadratureOptions = [&quadratureOrders](const std::string name,
          double relDist, int singleOrder, int doubleOrder) {


      ParameterList& paramList = quadratureOrders.sublist(name);
      paramList.set("maxRelDist",
              static_cast<double>(relDist),
              "(int) (Relative) distance of quadrature point to element.");

      paramList.set("singleOrder",
              static_cast<int>(singleOrder),
              "(int) (Relative) order of single regular integrals.");

      paramList.set("doubleOrder",
              static_cast<int>(doubleOrder),
              "(int) (Relative) order of double regular integrals.");

  };

  createQuadratureOptions("near",2,3,3);
  createQuadratureOptions("medium",4,2,2);
  createQuadratureOptions("far",std::numeric_limits<double>::infinity(),2,2);


  quadratureOrders.sublist("far").remove("maxRelDist");

  ParameterList& hmatParameters = parameters.sublist("HMat");

  hmatParameters.set("HMatAssemblyMode", std::string("GlobalAssembly"),
                     "(string) Specifies assembly mode. Allowed values are "
                     "GlobalAssembly and LocalAssembly");

  hmatParameters.set(
      "minBlockSize", static_cast<int>(50),
      "(int) Specifies the minimum block size below which blocks are "
      " assumed to be dense");

  hmatParameters.set(
      "maxBlockSize", static_cast<int>(2048),
      "(int) Specifies the maximum size of an admissible block.");

  hmatParameters.set("eta", static_cast<double>(1.2),
                     "(double) Specifies the block separation parameter eta");

  Teuchos::writeParameterListToXmlFile(parameters,"parameters.xml");
 
  return parameters;
}
}
