#Copyright (C) 2011 by the BEM++ Authors
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

# BEM++ setup script

import os, sys, traceback, platform
sys.path.append("installer")

from py_modules.tools import writeOptions, setDefaultConfigOption, pythonInfo, checkCreateDir, testBlas, testLapack, cleanUp, checkDeleteFile, checkInstallUpdates, installUpdates, normalizePath, to_bool, which
from ConfigParser import ConfigParser
from optparse import OptionParser




import py_modules.ahmed as ahmed
import py_modules.armadillo as armadillo
import py_modules.bempp as bempp
import py_modules.boost as boost
import py_modules.dune as dune
import py_modules.cmake as cmake
import py_modules.mkl as mkl
import py_modules.swig as swig
import py_modules.tbb as tbb
import py_modules.tools as tools
import py_modules.trilinos as trilinos

libraries = {'ahmed':ahmed,
             'armadillo':armadillo,
             'boost':boost,
             'dune':dune,
             'mkl':mkl,
             'swig':swig,
             'tbb':tbb,
             'trilinos':trilinos,
            }
library_names = sorted(libraries.keys())
library_names = sorted(library_names, key=lambda n: -1 if n == 'mkl' or n == 'swig' else 1)

###########################

# The following functions determine the directory in which the file resides.
# Trick taken from http://stackoverflow.com/questions/2632199/how-do-i-get-the-path-of-the-current-executed-file-in-python

def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def module_path():
    if we_are_frozen():
        return os.path.split(os.path.abspath(sys.executable))[0]
    return os.path.split(os.path.abspath(__file__))[0]

def downloadDependencies(root,config):

    print "Downloading dependencies"

    for dep in library_names:
        libraries[dep].download(root,config)

def createDirectories(root,config):

    dep_build_dir=config.get('Main','dependency_build_dir')
    checkCreateDir(dep_build_dir)
    dep_download_dir=config.get('Main','dependency_download_dir')
    checkCreateDir(dep_download_dir)
    bempp_build_dir=config.get('Bempp','build_dir')
    checkCreateDir(bempp_build_dir)
    checkCreateDir(root+"/installer/files")

    prefix = config.get('Main','prefix')
    checkCreateDir(prefix+"/bempp")
    checkCreateDir(prefix+"/bempp/lib")
    checkCreateDir(prefix+"/bempp/include")


def bootstrap(root,config):


    # Check for OS X Mavericks
    have_mavericks = False
    import platform
    plat = platform.system()
    if plat == 'Darwin':
        have_mavericks = True if platform.mac_ver()[0]=='10.9' else False

    cleanUp(root,config)
    createDirectories(root,config)

    prefix = config.get('Main','prefix')
    cmake_path = prefix+'/bempp/bin/cmake'
    if config.has_option('CMake','exe'):
	cmake_executable = tools.which(config.get('CMake','exe'))
	if cmake_executable is None:
	    raise Exception("CMake command specified in [CMake] section not found.")
        else:
            tools.setDefaultConfigOption(config,"CMake","exe",cmake_executable,overwrite=True)
    else:
        if not have_mavericks:
            print "CMake not found. Downloading and installing CMake..."
            cmake.download(root,config)
            cmake.prepare(root,config)
            cmake.configure(root,config)
            cmake.build(root,config)
            cmake.install(root,config)
            tools.setDefaultConfigOption(config,"CMake","exe",cmake_path,overwrite=True)
        else:
            raise Exception("On OS X Mavericks 'cmake' must be manually specified in the config file.") 

    downloadDependencies(root,config)

def prepareDependencies(root,config):


    for dep in library_names:
        libraries[dep].prepare(root,config)

    # Add soft link to Python library

    prefix = config.get('Main','prefix')
    py_lib = config.get('Python','lib')
    py_lib_name = os.path.basename(py_lib)
    py_soft_link_path = os.path.join(prefix,'bempp/lib/'+py_lib_name)
    if not os.path.lexists(py_soft_link_path):
	    os.symlink(py_lib,py_soft_link_path)


def installDependencies(root,config):

    for dep in library_names:
        libraries[dep].configure(root,config)
        libraries[dep].build(root,config)
        libraries[dep].install(root,config)


def prepare(root,config):
    # Test whether the main options are present
    if not config.has_option('Main','prefix'): raise Exception('prefix not defined')
    setDefaultConfigOption(config,'Main','cc','gcc')
    setDefaultConfigOption(config,'Main','cxx','g++')
    setDefaultConfigOption(config,'Main','architecture','intel64')
    setDefaultConfigOption(config,'Main','flags',"")
    setDefaultConfigOption(config,'Main','libs',"")
    setDefaultConfigOption(config,'Main','cflags',"")
    setDefaultConfigOption(config,'Main','cxxflags',"")
    setDefaultConfigOption(config,'Main','root_dir',root)
    setDefaultConfigOption(config,'Main','build_jobs',"1")

    # Retrieve path to configuration file
    optfile = config.get('Main','optfile')

    # Retrieve build directory
    setDefaultConfigOption(config,'Main','build_dir',root+'/build')
    build_dir = normalizePath(config, config.get('Main','build_dir'))
    # Set build directories for BEM++ and its dependencies
    config.set('Main','build_dir',build_dir)
    config.set('Bempp','build_dir',build_dir+'/bempp')
    config.set('Main','dependency_build_dir',build_dir+'/contrib')
    # Set
    config.set('Main','dependency_download_dir',root+'/installer/files')

    # Set default MKL/libs option
    setDefaultConfigOption(config,'MKL','lib',"-lmkl_rt")

    # Set default MPI options

    setDefaultConfigOption(config,'MPI','enable_mpi','false')
    setDefaultConfigOption(config,'MPI','mpi_cxx_libs','')
    setDefaultConfigOption(config,'MPI','mpi_include_dir','')
    enable_mpi = to_bool(config.get('MPI','enable_mpi'))

    if enable_mpi:
        mpi_include_dir = config.get('MPI','mpi_include_dir')
        cflags = config.get('Main','cflags')
        cxxflags = config.get('Main','cxxflags')
        config.set('Main','with_mpi','ON')
        if len(mpi_include_dir)>0:
	        config.set('Main','cflags',cflags+" -I"+mpi_include_dir)
	        config.set('Main','cxxflags',cxxflags+" -I"+mpi_include_dir)
    else:
        config.set('MPI','with_mpi','OFF')

 


    # Set empty BLAS/Lapack options if none exist
    setDefaultConfigOption(config,'BLAS','lib',"")
    setDefaultConfigOption(config,'LAPACK','lib',"")

    # Add the correct architecture parameters
    cflags = config.get('Main','cflags')
    cxxflags = config.get('Main','cxxflags')

    arch = config.get('Main','architecture')
    if not arch in ['ia32','ia64','intel64']:
        raise Exception("Architecture '"+arch+"' is not supported. "
                        "Supported architectures: ia32, ia64, intel64.")

    if sys.platform.startswith('darwin'):
        if arch=='intel64':
            param = '-arch x86_64'
        else:
            param = '-arch i386'
        config.set('Main','cflags',cflags+" "+param)
        config.set('Main','cxxflags',cxxflags+" "+param)
        setDefaultConfigOption(config,'Main','optflags','-O3 -march=core2')
    elif sys.platform.startswith('linux'):
        if arch=='intel64' or arch=='ia64':
            param = '-m64'
        else:
            param = '-m32'
        config.set('Main','cflags',cflags+" "+param)
        config.set('Main','cxxflags',cxxflags+" "+param)
        setDefaultConfigOption(config,'Main','optflags','-O3 -march=native')
    else:
        raise Exception("Platform '"+sys.platform+"' is not supported")

    # Add the correct Python options

    import numpy

    (py_exe,py_lib,py_include) = pythonInfo(config)
    setDefaultConfigOption(config,'Python','exe',py_exe)
    setDefaultConfigOption(config,'Python','lib',py_lib)
    setDefaultConfigOption(config,'Python','include_dir',py_include)
    setDefaultConfigOption(config,'Python','numpy_include_dir',numpy.get_include())


    # Check for OS X Mavericks
    have_mavericks = False
    plat = platform.system()
    if plat == 'Darwin':
        have_mavericks = True if platform.mac_ver()[0]=='10.9' else False
        # Add the CMake configuration

    if config.has_option('CMake','exe'):
        cmake_executable = tools.which(config.get('CMake','exe'))
        if cmake_executable is None:
            raise Exception("CMake command specified in [CMake] section not found.")
        tools.setDefaultConfigOption(config,'CMake','exe',cmake_executable,overwrite=True)
    else:
        # CMake must have been or will be  downloaded by the bootstrap mechanism
        cmake_executable = prefix+'/bempp/bin/cmake'
        if os.path.isfile(cmake_executable):
            tools.setDefaultConfigOption(config,'CMake','exe',cmake_executable,overwrite=True)

###########################

if __name__ == "__main__":

    usage = "usage: %prog [options] configuration_file"
    parser = OptionParser()
    parser.add_option("-b","--bootstrap", action="store_true",
                      help="Download dependencies and prepare directories")
    parser.add_option("-c", "--configure", action="store_true",
                      help="Configure the setup program")
    parser.add_option("-u","--update", action="store_true",help="Automatically update BEM++")
    parser.add_option("","--resume-update", action="store_true",help="Resume update process after changing the source tree")
    parser.add_option("-i", "--install", type="string", metavar="WHAT",
                      help="Build and install WHAT. Possible values for WHAT: "
                      "all (BEM++ and its dependencies), bempp (BEM++ only), " +
                      ", ".join(library_names) + " (particular BEM++ dependencies)")
    (options,args) = parser.parse_args()
    root=module_path()
    config=ConfigParser()
    if len(args) != 1:
        parser.error("Configuration file not specified")
    optfile = args[0]
    optfile_generated = root+"/"+os.path.basename(optfile)+".generated"
    try:
        optfileobj = open(optfile)
        config.readfp(optfileobj)
        optfileobj.close()
        optfile_full = os.path.abspath(os.path.expanduser(optfile))
        tools.setDefaultConfigOption(config,'Main','optfile',
                                     optfile_full,overwrite=True)
        prefix = normalizePath(config, config.get('Main','prefix'))
        tools.setDefaultConfigOption(config,'Main','prefix',prefix,
                                     overwrite=True)
    except Exception, e:
        print ("Parsing of configuration file '" + optfile +
               "' failed with error message:\n" + str(e))
        sys.exit(1)
    try:
        prepare(root,config)
        if options.resume_update:
            # Must be the first "if": the intention is that if this option
            # is present, the update procedure is resumed and all other work
            # modes are ignored.
            config = ConfigParser()
            if not os.path.exists(optfile_generated):
                print ("You must first successfully run bempp_setup.py "
                       "with the --configure (-c) option.")
                sys.exit(1)
            config.read(optfile_generated)
            installUpdates(root,config)
            sys.exit(0) # don't do anything else
        if options.update:
            config = ConfigParser()
            if not os.path.exists(optfile_generated):
                print ("You must first successfully run bempp_setup.py "
                       "with the --configure (-c) option.")
                sys.exit(1)
            config.read(optfile_generated)
            checkInstallUpdates(root,config)
        if options.bootstrap:
            bootstrap(root,config)
        if options.configure:
            checkDeleteFile(optfile_generated)
            prepareDependencies(root,config)
            bempp.prepare(root,config)
            testBlas(root,config)
            testLapack(root,config)
            opt_fp = open(optfile_generated,'w')
            config.write(opt_fp)
            opt_fp.close()
            print "Updated configuration written to '"+optfile_generated+"'"
            enable_mkl = tools.to_bool(config.get('MKL','enable_mkl'))
            if not enable_mkl:
                print ("----------------------------------------------------------\n"
                       "You configured BEM++ to use another BLAS and LAPACK\n"
                       "libraries than Intel MKL. For optimum performance, ensure\n"
                       "that your BLAS and LAPACK libraries are configured to work\n"
                       "in single-threaded mode, as otherwise threads spawned by\n"
                       "BLAS and LAPACK will compete for resources with those\n"
                       "spawned by BEM++. For instance, if you are using\n"
                       "GotoBLAS, set the environmental variable GOTO_NUM_THREADS\n"
                       "to '1' before running any programs using BEM++.\n")

        if options.install:
            config = ConfigParser()
            if not os.path.exists(optfile_generated):
                print ("You must first successfully run bempp_setup.py "
                       "with the --configure (-c) option.")
                sys.exit(1)
            config.read(optfile_generated)
            writeOptions(root,config)
            if options.install in library_names:
                libraries[options.install].configure(root,config)
                libraries[options.install].build(root,config)
                libraries[options.install].install(root,config)
            elif options.install == "all":
                installDependencies(root,config)
                bempp.configure(root,config)
                bempp.build(root,config)
                bempp.install(root,config)
            elif options.install == "bempp":
                bempp.configure(root,config)
                bempp.build(root,config)
                bempp.install(root,config)
            else:
                raise Exception("Library name not recognized.")
    except Exception, e:
        print "=" * 78, "\nBEM++ INSTALLATION FAILED WITH ERROR MESSAGE: \n"+ str(e)
        try:
            error_log = open(root+"/bempp_setup.err", "w")
            error_log.write(traceback.format_exc())
            error_log.close()
            print ("-" * 78 + "\n"
                   "Note: A stack trace of the BEM++ installer, providing detailed "
                   "information\non where the problem occured, has been written "
                   "to the file 'bempp_setup.err'.\nPlease include this "
                   "file if you report the problem to BEM++ developers.")
        except:
           pass
        sys.exit(1)
