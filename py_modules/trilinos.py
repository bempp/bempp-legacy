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


trilinos_fname='trilinos-10.10.2-Source.tar.gz'
trilinos_extract_dir='trilinos-10.10.2-Source'
trilinos_dir='trilinos'
trilinos_url='http://trilinos.sandia.gov/download/files/trilinos-10.10.2-Source.tar.gz'

def download(root,config):
    tools.download(trilinos_fname,trilinos_url,root+"/contrib/files")


def prepare(root,config):

    trilinos_full_dir=root+"/contrib/"+trilinos_dir
    trilinos_download_name=root+"/contrib/files/"+trilinos_fname

    prefix=config.get('Main','prefix')

    tools.checkDeleteDirectory(trilinos_full_dir)
    
    print "Extracting Trilinos"
    tools.extract_file(root+"/contrib/files/"+trilinos_fname,root+"/contrib/")
    os.rename(root+"/contrib/"+trilinos_extract_dir,root+"/contrib/"+trilinos_dir)
    shutil.copy(root+"/contrib/build_scripts/posix/trilinos_build.sh",trilinos_full_dir+"/trilinos_build.sh")
    print "Patching ..."
    patch=py_patch.fromfile(root+"/contrib/patch/Thyra_BelosLinearOpWithSolve_def.patch")
    cwd=os.getcwd()
    os.chdir(root+"/contrib/trilinos//packages/stratimikos/adapters/belos/src")
    patch.apply()
    os.chdir(cwd)
        
    tools.setDefaultConfigOption(config,'Trilinos','cmake_path',prefix+"/bempp/lib/cmake/Trilinos/",overwrite=True)

def configure(root,config):

    trilinos_full_dir=root+"/contrib/"+trilinos_dir
    trilinos_download_name=root+"/contrib/files/"+trilinos_fname

    print "Configuring Trilinos"
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir)
    tools.checkDeleteDirectory(trilinos_full_dir+"/build")
    subprocess.check_call("sh ./trilinos_build.sh",shell=True)
    os.chdir(cwd)

def build(root,config):

    trilinos_full_dir=root+"/contrib/"+trilinos_dir

    print "Build Trilinos"
    njobs = tools.to_int(config.get('Main','build_jobs'))
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir+"/build")
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

def install(root,config):
    print "Install Trilinos"
    trilinos_full_dir=root+"/contrib/"+trilinos_dir
    cwd=os.getcwd()
    os.chdir(trilinos_full_dir+"/build")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)

    
    
        
        
        
            

        
