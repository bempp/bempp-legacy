%{
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
  %}


%warnfilter(401) Bempp::StandardLocalAssemblerFactoryForOperatorsOnSurfaces;

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

class StandardLocalAssemblerFactoryForOperatorsOnSurfaces(Template2,ScalarSpace):
    """Standard assembler factory"""
    def __init__(self,dtype1,dtype2,*args,**kwargs):
          super(StandardLocalAssemblerFactoryForOperatorsOnSurfaces,self).__init__('StandardLocalAssemblerFactoryForOperatorsOnSurfaces',dtype1,dtype2,*args,**kwargs)

	  %}
