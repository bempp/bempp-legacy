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
from py_modules import python_patch as py_patch


trilinos_fname='trilinos-11.2.3-Source.tar.gz'
trilinos_extract_dir='trilinos-11.2.3-Source'
trilinos_dir='trilinos'
trilinos_url='http://trilinos.sandia.gov/download/files/trilinos-11.2.3-Source.tar.gz'

def download(root,config,force=False):
    dep_download_dir=config.get('Main','dependency_download_dir')
    tools.download(trilinos_fname,trilinos_url,dep_download_dir,force)

def prepare(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    dep_download_dir=config.get('Main','dependency_download_dir')

    trilinos_full_dir=dep_build_dir+"/"+trilinos_dir
    trilinos_download_name=dep_download_dir+"/"+trilinos_fname

    prefix=config.get('Main','prefix')

    tools.checkDeleteDirectory(trilinos_full_dir)

    print "Extracting Trilinos"
    try:
        tools.extract_file(dep_download_dir+"/"+trilinos_fname,dep_build_dir)
    except IOError:
        # Possibly a corrupted/truncated file. Try to download once again
        download(root,config,force=True)
        tools.extract_file(dep_download_dir+"/"+trilinos_fname,dep_build_dir)
    os.rename(dep_build_dir+"/"+trilinos_extract_dir,
              dep_build_dir+"/"+trilinos_dir)
    shutil.copy(root+"/installer/build_scripts/posix/trilinos_build.sh",
                trilinos_full_dir+"/trilinos_build.sh")
    print "Patching Trilinos"
    cwd=os.getcwd()
    os.chdir(dep_build_dir+"/trilinos/packages/stratimikos/adapters/belos/src")
    patch=py_patch.fromfile(root+"/installer/patches/Thyra_BelosLinearOpWithSolve_def.patch")
    patch.apply()
    os.chdir(dep_build_dir+"/trilinos/packages/thyra/core/src/support/nonlinear/model_evaluator/client_support")
    patch=py_patch.fromfile(root+"/installer/patches/thyra_static_initialization_order.patch")
    patch.apply()
    os.chdir(dep_build_dir+"/trilinos/cmake/tribits/modules")
    patch=py_patch.fromfile(root+"/installer/patches/trilinos_find_python_interp.patch")
    patch.apply()
    os.chdir(dep_build_dir+"/trilinos/packages/teuchos/numerics/src")
    patch=py_patch.fromfile(root+"/installer/patches/Teuchos_LAPACK.hpp.patch")
    patch.apply()
    patch=py_patch.fromfile(root+"/installer/patches/Teuchos_LAPACK.cpp.patch")
    patch.apply()
    os.chdir(dep_build_dir+"/trilinos/packages/epetra/src")
    patch=py_patch.fromfile(root+"/installer/patches/Epetra_ConfigDefs.h.patch")
    patch.apply()
    os.chdir(cwd)

    tools.setDefaultConfigOption(config,'Trilinos','cmake_path',prefix+"/bempp/lib/cmake/Trilinos/",overwrite=True)

    tools.setCompilerOptions(config,'Trilinos')
   

def configure(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    trilinos_full_dir=dep_build_dir+"/"+trilinos_dir

    print "Configuring Trilinos"
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir)
    tools.checkDeleteDirectory(trilinos_full_dir+"/build")
    subprocess.check_call("sh ./trilinos_build.sh",shell=True)
    os.chdir(cwd)

def build(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    trilinos_full_dir=dep_build_dir+"/"+trilinos_dir

    print "Build Trilinos"
    njobs = tools.to_int(config.get('Main','build_jobs'))
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir+"/build")
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

def install(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    trilinos_full_dir=dep_build_dir+"/"+trilinos_dir
    print "Install Trilinos"
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir+"/build")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)
