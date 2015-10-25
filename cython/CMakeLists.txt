# Figures out where to install python stuff

add_custom_target(cython-headers)
add_custom_target(cython)
add_subdirectory(bempp)
add_dependencies(cython cython-headers)

