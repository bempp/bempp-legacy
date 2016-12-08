# Simple function to call python
function(call_python OUTPUT)
   # First tries adding prefix
   execute_process(
     COMMAND ${PYTHON_EXECUTABLE} -c "${ARGN}"
     RESULT_VARIABLE result
     OUTPUT_VARIABLE output
   )
   if(result EQUAL 0)
       string(STRIP "${output}" output)
       set(${OUTPUT} "${output}" PARENT_SCOPE)
   endif()
endfunction()
