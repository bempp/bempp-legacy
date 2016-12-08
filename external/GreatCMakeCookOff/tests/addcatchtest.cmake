option(tests          "Enable testing."                         on)

if(tests)
  list(APPEND CMAKE_MODULE_PATH ${cookoff_path}/scripts ${cookoff_path}/modules)
  include(AddCatchTest)
  enable_testing()
endif(tests)

file(WRITE ${CMAKE_SOURCE_DIR}/mytest.cc
     "#include <catch.hpp>\n\n"
     "TEST_CASE(\"A test case\") {\n"
     "  CHECK(true);\n"
     "}\n\n"
)


add_catch_test(mytest)
