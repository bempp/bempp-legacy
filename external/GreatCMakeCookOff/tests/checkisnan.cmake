list(APPEND CMAKE_MODULE_PATH ${cookoff_path}/scripts ${cookoff_path}/modules)
include(CheckIsNaN)
if(NOT ISNAN_VARIATION)
  message(STATUS "Could not find working isnan.")
endif(NOT ISNAN_VARIATION)

enable_language(CXX)
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/isnan.cc.in
     "\@ISNAN_HEADERS\@\n"
     "#include <limits>\n"
     "#include <exception>\n"
     "#define not_a_number(X) \@ISNAN_VARIATION\@(X)\n"
     "int main() {\n"
     "  if(not_a_number(0.5e0)) throw std::exception(); \n"
     "  if(not_a_number(2l)) throw std::exception(); \n"
     "  if(not not_a_number(std::numeric_limits<double>::quiet_NaN())) throw std::exception(); \n"
     "  return 0;\n"
     "}\n")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/isnan.cc.in ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
