enable_language(C)
function(check_macro_existence NAME OUTVAR)
    set(bindir "${PROJECT_BINARY_DIR}/CMakeFiles/check_macro_existence")
    file(WRITE "${bindir}/${NAME}_positive.c"
        "#ifdef ${NAME}\n"
        "# error THIS IS AN ERROR\n"
        "#endif\n"
        "int main(int argc, char const *argv[]) { return 0; }\n"
    )
    file(WRITE "${bindir}/${NAME}_negative.c"
        "#ifndef ${NAME}\n"
        "# error THIS IS AN ERROR\n"
        "#endif\n"
        "int main(int argc, char const *argv[]) { return 0; }\n"
    )
    try_compile(POSITIVE "${bindir}" "${bindir}/${NAME}_positive.c")
    try_compile(NEGATIVE "${bindir}" "${bindir}/${NAME}_negative.c")
    if(NOT POSITIVE AND NEGATIVE)
        set(${OUTVAR} TRUE PARENT_SCOPE)
    elseif(POSITIVE AND NOT NEGATIVE)
        set(${OUTVAR} FALSE PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Could not determine existence of ${NAME}.")
    endif()
endfunction()

if(NOT DETECTED_ARCHITECTURE AND APPLE AND CMAKE_OSX_ARCHITECTURES)
    set(DETECTED_ARCHITECTURE "${CMAKE_OSX_ARCHITECTURES}" CACHE STRING
        "Type of the target architecture")
endif()
if(NOT DETECTED_ARCHITECTURE)
    foreach(macroname __ia64 __ia64__ _M_IA64)
        check_macro_existence(${macroname} DOES_EXIST)
        if(DOES_EXIST)
            set(DETECTED_ARCHITECTURE "ia64" CACHE STRING
                "Type of the target architecture")
            break()
        endif()
    endforeach()
endif()
if(NOT DETECTED_ARCHITECTURE)
    foreach(macroname __x86_64 __x86_64__ __amd64 _M_X64)
        check_macro_existence(${macroname} DOES_EXIST)
        if(DOES_EXIST)
            set(DETECTED_ARCHITECTURE "x86_64" CACHE STRING
                "Type of the target architecture")
            break()
        endif()
    endforeach()
endif()
if(NOT DETECTED_ARCHITECTURE)
    foreach(macroname __i386 __i386__ _M_IX86)
        check_macro_existence(${macroname} DOES_EXIST)
        if(DOES_EXIST)
            set(DETECTED_ARCHITECTURE "i386" CACHE STRING
                "Type of the target architecture")
            break()
        endif()
    endforeach()
endif()
