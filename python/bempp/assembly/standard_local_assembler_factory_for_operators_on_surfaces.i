%{
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
  %}


%warnfilter(401) Bempp::StandardLocalAssemblerFactoryForOperatorsOnSurfaces;

%feature("pythonappend") Bempp::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<float,float>::StandardLocalAssemblerFactoryForOperatorsOnSurfaces(const AccuracyOptions& accuracyOptions)%{
self._BasisFunctionType='float32'
self._ResultType='float32'
  %}





%include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"


namespace Bempp
{

  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_float32_float32) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<float,float>;
  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_float32_complex64) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<float,std::complex<float> >;
  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_float64_float64) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double,double>;
  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_float64_complex128) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double,std::complex<double> >;
  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_complex64_complex64) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<float>,std::complex<float> >;
  %template(StandardLocalAssemblerFactoryForOperatorsOnSurfaces_complex128_complex128) StandardLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<double>,std::complex<double> >;
 
    }



%pythoncode %{

  class StandardLocalAssemblerFactoryForOperatorsOnSurfaces(BasisFunctionType='float64',ResultType='float64',AccuracyOptions=None):
      if ResultType is not None:
            dtype=checkType(

	  %}
