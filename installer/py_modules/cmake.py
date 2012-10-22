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

import os,urllib,shutil,subprocess,sys
from py_modules import tools


cmake_fname='cmake-2.8.9.tar.gz'
cmake_url='http://www.cmake.org/files/v2.8/cmake-2.8.9.tar.gz'
cmake_extract_dir='cmake-2.8.9'
cmake_dir='cmake'

def download(root,config,force=False):
    dep_download_dir=config.get('Main','dependency_download_dir')
    cmake_download_name=dep_download_dir+"/"+cmake_fname
    tools.download(cmake_fname,cmake_url,dep_download_dir,force)

def prepare(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    dep_download_dir=config.get('Main','dependency_download_dir')
    prefix=config.get('Main','prefix')

    cmake_full_dir=dep_build_dir+"/"+cmake_dir
    cmake_download_name=dep_download_dir+"/"+cmake_fname

    cmake_executable=prefix+"/bempp/bin/cmake"

    tools.checkDeleteDirectory(cmake_full_dir)

    print "Extracting CMake"
    try:
        tools.extract_file(dep_download_dir+"/"+cmake_fname,dep_build_dir)
    except IOError:
        # Possibly a corrupted/truncated file. Try to download once again
        download(root,config,force=True)
        tools.extract_file(dep_download_dir+"/"+cmake_fname,dep_build_dir)
    os.rename(dep_build_dir+"/"+cmake_extract_dir,cmake_full_dir)


def configure(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    cmake_full_dir=dep_build_dir+"/"+cmake_dir
    prefix=config.get('Main','prefix')
    cmake_prefix=prefix+"/bempp"
    print "Configuring CMake"
    cwd=os.getcwd()
    os.chdir(cmake_full_dir)
    njobs = tools.to_int(config.get('Main','build_jobs'))
    subprocess.check_call(["./bootstrap", "--prefix="+cmake_prefix,"--parallel="+str(njobs)])
    os.chdir(cwd)

def build(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    cmake_full_dir=dep_build_dir+"/"+cmake_dir
    njobs = tools.to_int(config.get('Main','build_jobs'))
    print "Building CMake"
    cwd=os.getcwd()
    os.chdir(cmake_full_dir)
    subprocess.check_call(["make", "-j"+str(njobs)])
    os.chdir(cwd)

def install(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    cmake_full_dir=dep_build_dir+"/"+cmake_dir
    print "Installing CMake"
    cwd=os.getcwd()
    os.chdir(cmake_full_dir)
    subprocess.check_call(["make", "install"])
    os.chdir(cwd)
