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

# This script downloads the third-party libraries required by BEM++


import sys,os
from py_modules.tools import writeOptions, setDefaultConfigOption, pythonInfo, checkCreateDir
from ConfigParser import ConfigParser
from optparse import OptionParser


import py_modules.boost as boost
import py_modules.armadillo as armadillo
import py_modules.tbb as tbb
import py_modules.dune as dune
import py_modules.trilinos as trilinos
import py_modules.bempp as bempp
import py_modules.ahmed as ahmed
import py_modules.mkl as mkl

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

def downloadAll(root,config):

    checkCreateDir(root+"/contrib/files")

    armadillo.download(root,config)
    tbb.download(root,config)
    mkl.download(root,config)
    ahmed.download(root,config)
    boost.download(root,config)
    trilinos.download(root,config)
    dune.download(root,config)
    

def prepareAll(root,config):

    prefix=config.get('Main','prefix')
    checkCreateDir(prefix+"/bempp")
    checkCreateDir(prefix+"/bempp/lib")
    checkCreateDir(prefix+"/bempp/include")

    #armadillo.prepare(root,config)
    #tbb.prepare(root,config)
    mkl.prepare(root,config)
    ahmed.prepare(root,config)
    #boost.prepare(root,config)
    #trilinos.prepare(root,config)
    #dune.prepare(root,config)
    
def configureAll(root,config):
    
    mkl.configure(root,config)
    #boost.configure(root,config)
    #dune.configure(root,config)
    #tbb.configure(root,config)
    #armadillo.configure(root,config)
    #trilinos.configure(root,config)
    ahmed.configure(root,config)
    #bempp.configure(root,config)
    pass

def buildAll(root,config):

    mkl.build(root,config)
    #boost.build(root,config)
    #tbb.build(root,config)
    #armadillo.build(root,config)
    #trilinos.build(root,config)
    ahmed.build(root,config)
    #dune.build(root,config)
    #bempp.build(root,config)
    pass

def installAll(root,config):

    mkl.install(root,config)
    #dune.install(root,config)
    #tbb.install(root,config)
    #armadillo.install(root,config)
    #boost.install(root,config)
    #trilinos.install(root,config)
    ahmed.install(root,config)
    

def prepare(root,config):
    # Test whether the main options are present
    if not config.has_option('Main','prefix'): raise Exception('prefix not defined')
    setDefaultConfigOption(config,'Main','cc','gcc')
    setDefaultConfigOption(config,'Main','cxx','g++')
    setDefaultConfigOption(config,'Main','architecture','intel64')
    setDefaultConfigOption(config,'Main','cflags',"")
    setDefaultConfigOption(config,'Main','cxxflags',"")
    setDefaultConfigOption(config,'Main','root_dir',root)
    setDefaultConfigOption(config,'Main','build_jobs',1)

    # Add the correct architecture parameters
    cflags = config.get('Main','cflags')
    cxxflags = config.get('Main','cxxflags')

    arch = config.get('Main','architecture')
    if not arch in ['intel64','i386']: raise Exception('Architecture not supported.')
    
    if sys.platform.startswith('darwin'):
        if arch=='intel64':
            param = '-arch x86_64'
        else:
            param = '-arch i386'
        config.set('Main','cflags',cflags+" "+param)
        config.set('Main','cxxflags',cxxflags+" "+param)
    elif sys.platform.startswith('linux'):
        if arch=='x64':
            param = '-m64'
        else:
            param = '-m32'
        config.set('Main','cflags',cflags+" "+param)       
        config.set('Main','cxxflags',cxxflags+" "+param)
    else:
        raise Exception("Platform not supported")

    # Add the correct Python options

    (py_exe,py_lib,py_include) = pythonInfo()
    setDefaultConfigOption(config,'Python','exe',py_exe)
    setDefaultConfigOption(config,'Python','lib',py_lib)
    setDefaultConfigOption(config,'Python','include_dir',py_include)
     
###########################

if __name__ == "__main__":
 
    parser = OptionParser()
    parser.add_option("-c", "--configure", action="store_true", dest="configure", default=False)
    parser.add_option("-b", "--build", action="store_true", dest="build", default=False)
    parser.add_option("-d", "--download", action="store_true", dest="download", default=False)
    parser.add_option("-p", "--prepare", action="store_true", dest="prepare", default=False)
    parser.add_option("-i", "--install", action="store_true", dest="install", default=False)
    (options,args) = parser.parse_args()
    root=module_path()
    config=ConfigParser()
    optfile = args[0]
    optfile_updated = optfile+".new"
    config.read(optfile)
    prepare(root,config)
    if options.download:
        downloadAll(root,config)
    if options.prepare:
        prepareAll(root,config)
        opt_fp = open(optfile_updated,'w')
        config.write(opt_fp)
        opt_fp.close()
        writeOptions(root,config)
    if options.configure:
        config = ConfigParser()
        config.read(optfile_updated)
        configureAll(root,config)
    if options.build:
        config = ConfigParser()
        config.read(optfile_updated)
        buildAll(root,config)
    if options.install:
        config = ConfigParser()
        config.read(optfile_updated)
        installAll(root,config)


    
