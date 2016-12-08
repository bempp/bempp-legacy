# Changes/unchanges the root path
# Makes it possible to look for packages stricly inside the root
macro(_set_root_path PATH TYPE)
    foreach(save MODE_PACKAGE MODE_INCLUDE MODE_PROGRAM MODE_LIBRARY)
        set(_save_root_path_${save} ${CMAKE_FIND_ROOT_PATH_${save}})
        set(CMAKE_FIND_ROOT_PATH_${save} ${TYPE})
    endforeach()
    set(_save_root_path "${CMAKE_FIND_ROOT_PATH}")
    set(CMAKE_FIND_ROOT_PATH "${PATH}")
endmacro()

# Unchanges the root path
macro(_unset_root_path)
    foreach(save MODE_PACKAGE MODE_INCLUDE MODE_PROGRAM MODE_LIBRARY)
        set(CMAKE_FIND_ROOT_PATH_${save} ${_save_root_path_${save}})
    endforeach()
    set(CMAKE_FIND_ROOT_PATH "${_save_root_path}")
endmacro()


