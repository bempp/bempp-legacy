# Copyright (c) 2015, Damien Buhl 
# All rights reserved.
include(hunter_add_version)
include(hunter_download)
include(hunter_pick_scheme)
include(hunter_add_package)
include(hunter_configuration_types)

# Makes it possible to use syste cfitsio
hunter_add_version(
    PACKAGE_NAME
    GMP
    VERSION
    "6.1.1"
    URL
    "https://gmplib.org/download/gmp/gmp-6.1.1.tar.xz"
    SHA1
    4da491d63ef850a7662f41da27ad1ba99c2dbaa1
)

hunter_pick_scheme(DEFAULT GMP)
hunter_configuration_types(GMP CONFIGURATION_TYPES Release)
hunter_download(PACKAGE_NAME GMP)
