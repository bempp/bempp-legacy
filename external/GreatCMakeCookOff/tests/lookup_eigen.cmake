find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

include(PackageLookup)
lookup_package(Eigen3 DOWNLOAD_BY_DEFAULT KEEP)

if(EIGEN3_INCLUDE_DIR)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()
file(WRITE "${CMAKE_SOURCE_DIR}/main.cc"
    "#include <iostream>\n"
    "#include <Eigen/Dense>\n"
    "using Eigen::MatrixXd;\n"
    "int main() {\n"
    "  MatrixXd m(2,2);\n"
    "  m(0,0) = 3;\n"
    "  m(1,0) = 2.5;\n"
    "  m(0,1) = -1;\n"
    "  m(1,1) = m(1,0) + m(0,1);\n"
    "  MatrixXd msquared = m * m;\n"
    "  return 0;\n"
    "}\n"
)
add_executable(eigen "${CMAKE_SOURCE_DIR}/main.cc")
add_dependencies(eigen lookup_dependencies)
