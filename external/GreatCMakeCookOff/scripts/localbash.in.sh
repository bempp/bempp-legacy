#!@ENV_EXECUTABLE@ @BASH_EXECUTABLE@
if [ -e @PROJECT_BINARY_DIR@/paths/ldpaths ]; then
    add_to_ld() {
        if [ -d "$1" ] && [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
            LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+"$LD_LIBRARY_PATH:"}$1"
        fi
        if [ -d "$1" ] && [[ ":$DYLD_FALLBACK_LIBRARY_PATH:" != *":$1:"* ]]; then
            DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH:+"$DYLD_FALLBACK_LIBRARY_PATH:"}$1"
        fi
    }
    while read -r line; do
       add_to_ld $line
    done < @PROJECT_BINARY_DIR@/paths/ldpaths
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH
fi
if [ -e @PROJECT_BINARY_DIR@/paths/pypaths.pth ] && [ -n "@env_PYTHON@" ]; then
    add_to_py() {
      if [ -d "$1" ] && [[ ":$PYTHONPATH:" != *":$1:"* ]]; then
        PYTHONPATH="${PYTHONPATH:+"$PYTHONPATH:"}$1"
      fi
    }
    while read -r line; do
       add_to_py $line
    done < @PROJECT_BINARY_DIR@/paths/pypaths.pth
    export PYTHONPATH=$PYTHONPATH
fi
export PATH=@PROJECT_BINARY_DIR@/external/bin:$PATH
if [ -n "@env_WORKING_DIRECTORY@" ]; then
    cd @env_WORKING_DIRECTORY@
fi
if [ -n "@env_EXECUTABLE@" ]; then
    @env_EXECUTABLE@ "$@"
else
    eval $@
fi
