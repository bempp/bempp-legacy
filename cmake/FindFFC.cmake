# Find FFC

include(PythonPackageLookup)
include(CMakeParseArguments)
include(PythonPackage)
include(PackageLookup)


find_python_package(ffc)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFC
    "FFC could not be found. Downloading and installing FFC."
    ffc_FOUND)
