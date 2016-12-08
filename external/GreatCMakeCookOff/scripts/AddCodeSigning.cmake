# Adds the ability to do code signing for macos
# It avoids getting a pop-up all the time when running MPI codes.
# The process is as follows:
#
# Security nerds, avert your eyes.
#
# 1. create a *code-signing* certificate:
#     - open key-chain access
#     - select the login keychain (top-left panel)
#     - menu KeyChain Access/Certificate Assistant/Create a certificate
#     - choose a certificate name, self-signed, code-signing
#     - curse Apple for not providing a command-line utility
# 1. run cmake with:
#
#     ~~~bash
#     cmake -DCERTIFICATE="certificate name" .
#     ~~~
#
# 1. to every mpi target, add in the relevant CMakeLists.txt:
#
#     ~~~CMake
#     add_codesigning(targetname)
#     ~~~
#
# 1. the very first time the certificate is used (during build) click always allow

set(CERTIFICATE "" CACHE STRING "Name of the certificate to use for code-signing")
function(add_codesign target)
  if(APPLE AND CERTIFICATE)
    add_custom_command(
      TARGET ${target}
      POST_BUILD
      COMMAND codesign -f -s "${CERTIFICATE}" $<TARGET_FILE:${target}>
      COMMENT "Codsigning ${target}"
      VERBATIM
    )
  endif()
endfunction()
