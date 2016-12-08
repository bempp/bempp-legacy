# Tries and find which aligned alloc to use
# defines ALLOC_HEADER, which should be included in a configuration files somewher as #
#  @ALLOC_HEADER@, and ALLOC_VARIATION which is the name of the actual alloc.
#
# The simplest is to include in some header file something like
# 
# @ALIGNED_ALLOC_HEADER@
# inline void* my_aligned_alloc(size_t _size, size_t _alignment) { @ALIGNED_ALLOC_VARIATION@; }
# inline void my_aligned_free(void *_pointer) { @ALIGNED_FREE_VARIATION@; }
#
# The name of the arguments is not optional.

if(FOUND_ALIGNED_ALLOC)
  return()
endif()

include(CheckCXXSourceRuns)

macro(check_aligned_alloc NAME RESULT HEADER ALLOC FREE)
 CHECK_CXX_SOURCE_RUNS(
    "${HEADER}\n 
     inline void* my_aligned_alloc(size_t _size, size_t _alignment) { ${ALLOC}; }
     inline void my_aligned_free(void *_pointer) { ${FREE}; }\n
     int main() {
       void *ptr = my_aligned_alloc(10, 16);
       my_aligned_free(ptr);
       return 0;
     }\n" 
  aligned_alloc_${NAME})
  set(${RESULT} ${aligned_alloc_${NAME}})
endmacro()

set(ALIGNED_ALLOC_HEADER "#include <malloc.h>")
set(ALIGNED_ALLOC_VARIATION "return _aligned_malloc(_size, _alignment);")
set(ALIGNED_FREE_VARIATION "_aligned_free(_pointer);")
set(ALIGNED_ALLOC_NAME Windows)
check_aligned_alloc( Windows FOUND_ALIGNED_ALLOC
                     "${ALIGNED_ALLOC_HEADER}"
                     "${ALIGNED_ALLOC_VARIATION}"
                     "${ALIGNED_FREE_VARIATION}" )
if(NOT FOUND_ALIGNED_ALLOC)
  set(ALIGNED_ALLOC_HEADER "#include <stdlib.h>")
  set(ALIGNED_ALLOC_VARIATION
    "size_t const residual = _size % _alignment;"
    "size_t const actual_size = residual > 0 ? _size + _alignment - residual: _size;"
    "return ::aligned_alloc(_alignment, actual_size);")
  set(ALIGNED_FREE_VARIATION "free(_pointer);")
  set(ALIGNED_ALLOC_NAME C11)
  check_aligned_alloc( C11 FOUND_ALIGNED_ALLOC
                       "${ALIGNED_ALLOC_HEADER}"
                       "${ALIGNED_ALLOC_VARIATION}"
                       "${ALIGNED_FREE_VARIATION}" )
endif()
if(NOT FOUND_ALIGNED_ALLOC)
  set(ALIGNED_ALLOC_HEADER "#include <stdlib.h>")
  set(ALIGNED_ALLOC_VARIATION
    "size_t const residual = _size % _alignment;"
    "size_t const actual_size = residual > 0 ? _size + _alignment - residual: _size;"
    "void *result;"
    "if(posix_memalign(&result, _alignment, actual_size)) return NULL;"
    "return result;")
  set(ALIGNED_FREE_VARIATION "free(_pointer);")
  set(ALIGNED_ALLOC_NAME POSIX)
  check_aligned_alloc( POSIX FOUND_ALIGNED_ALLOC
                       "${ALIGNED_ALLOC_HEADER}"
                       "${ALIGNED_ALLOC_VARIATION}"
                       "${ALIGNED_FREE_VARIATION}" )
endif()
if(NOT FOUND_ALIGNED_ALLOC)

  set(ALIGNED_ALLOC_HEADER "#include <stdlib.h>")
  set(ALIGNED_ALLOC_VARIATION "size_t const residual = (_size + sizeof(void**)) % _alignment;"
      "size_t actual_size = _size + sizeof(void**);"
      "if(residual > 0) actual_size += _alignment - residual;"
      "void *const original = static_cast<void*>(new char[actual_size]);"
      "size_t const a = size_t(static_cast<void**>(original)+1) % _alignment;"
      "char * result = reinterpret_cast<char*>(static_cast<void**>(original)+1);"
      "if(a > 0) result += _alignment - a;"
      "*(reinterpret_cast<void**>(result) - 1) = original;"
      "return static_cast<void*>(result)")
  set(ALIGNED_ALLOC_NAME GreatCMakeCookProvides)
  set(ALIGNED_FREE_VARIATION "free(*(static_cast<void**>(_pointer)-1));")
  check_aligned_alloc( HomeMade FOUND_ALIGNED_ALLOC
                       "${ALIGNED_ALLOC_HEADER}"
                       "${ALIGNED_ALLOC_VARIATION}"
                       "${ALIGNED_FREE_VARIATION}" )
endif()

if(FOUND_ALIGNED_ALLOC)
  set(FOUND_ALIGNED_ALLOC ${FOUND_ALIGNED_ALLOC} CACHE BOOL
      "Found an aligned allocation function")
  set(ALIGNED_ALLOC_HEADER "${ALIGNED_ALLOC_HEADER}" CACHE STRING
      "Aligned allocation header")
  set(ALIGNED_ALLOC_VARIATION "${ALIGNED_ALLOC_VARIATION}" CACHE STRING
      "Aligned allocation function")
  set(ALIGNED_FREE_VARIATION "${ALIGNED_FREE_VARIATION}" CACHE STRING
      "Aligned free function")
  message(STATUS "[Aligned Alloc] Variation found in ${ALIGNED_ALLOC_NAME}")
else()
  set(FOUND_ALIGNED_ALLOC FALSE)
  set(ALIGNED_ALLOC_HEADER)
  set(ALIGNED_ALLOC_VARIATION)
  set(ALIGNED_FREE_VARIATION)
endif(FOUND_ALIGNED_ALLOC)
