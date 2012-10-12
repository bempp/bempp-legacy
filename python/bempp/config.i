%{
#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_alugrid.hpp"
#include "bempp/common/config_data_types.hpp"
#include "bempp/common/config_opencl.hpp"
#include "bempp/common/config_trilinos.hpp"
%}

%include "bempp/common/config_ahmed.hpp"
%include "bempp/common/config_alugrid.hpp"
%include "bempp/common/config_data_types.hpp"
%include "bempp/common/config_opencl.hpp"
%include "bempp/common/config_trilinos.hpp"

%inline %{
#ifdef WITH_AHMED
const bool _withAhmed = true;
#else
const bool _withAhmed = false;
#endif // WITH_AHMED
%}
