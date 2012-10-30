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


boost_fname='boost_1_49.tar.gz'
boost_url='http://www.bempp.org/files/boost_1_49.tar.gz'
boost_extract_dir='boost-cmake'
boost_dir='boost'
boost_version="1.49.0"

def download(root,config,force=False):
    dep_download_dir=config.get('Main','dependency_download_dir')
    tools.download(boost_fname,boost_url,dep_download_dir,force)

def prepare(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    boost_full_dir=dep_build_dir+"/"+boost_dir
    dep_download_dir=config.get('Main','dependency_download_dir')
    boost_download_name=dep_download_dir+"/"+boost_fname

    prefix=config.get('Main','prefix')
    boost_include_dir=prefix+"/bempp/include"

    if sys.platform.startswith('darwin'):
        unit_test_lib_name="libboost_unit_test_framework-mt.dylib"
    elif sys.platform.startswith('linux'):
        unit_test_lib_name="libboost_unit_test_framework-mt.so"
    else:
        raise Exception("Platform not supported")

    boost_unit_test_lib=prefix+"/bempp/lib/"+unit_test_lib_name

    tools.checkDeleteDirectory(boost_full_dir)

    print "Extracting Boost"
    try:
        tools.extract_file(boost_download_name,dep_build_dir)
    except IOError:
        # Possibly a corrupted/truncated file. Try to download once again
        download(root,config,force=True)
        tools.extract_file(boost_download_name,dep_build_dir)
    os.rename(dep_build_dir+"/"+boost_extract_dir,boost_full_dir)
    shutil.copy(root+"/installer/build_scripts/posix/boost_build.sh",
                boost_full_dir+"/boost_build.sh")

    tools.setDefaultConfigOption(config,"Boost","unit_test_lib",
                                 boost_unit_test_lib,overwrite=True)
    tools.setDefaultConfigOption(config,"Boost","include_dir",
                                 boost_include_dir,overwrite=True)

    tools.setCompilerOptions(config,'Boost')

def configure(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    boost_full_dir=dep_build_dir+"/"+boost_dir
    print "Configuring Boost"
    cwd=os.getcwd()
    os.chdir(boost_full_dir)
    tools.checkDeleteDirectory(boost_full_dir+"/build")
    subprocess.check_call("sh ./boost_build.sh",shell=True)
    os.chdir(cwd)


def build(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    boost_full_dir=dep_build_dir+"/"+boost_dir
    njobs=tools.to_int(config.get('Main','build_jobs'))
    print "Build Boost"
    cwd=os.getcwd()
    os.chdir(boost_full_dir+"/build")
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

def install(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    boost_full_dir=dep_build_dir+"/"+boost_dir
    prefix=config.get('Main','prefix')
    print "Install Boost"
    cwd=os.getcwd()
    os.chdir(boost_full_dir+"/build")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)









