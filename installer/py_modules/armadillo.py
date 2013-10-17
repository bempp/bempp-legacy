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


#arma_fname='armadillo-3.2.0.tar.gz'
#arma_url='http://sourceforge.net/projects/arma/files/armadillo-3.2.0.tar.gz'
#arma_extract_dir='armadillo-3.2.0'
arma_fname='armadillo-3.920.1.tar.gz'
arma_url='http://sourceforge.net/projects/arma/files/armadillo-3.920.1.tar.gz'
arma_extract_dir='armadillo-3.920.1'
arma_dir='armadillo'


def download(root,config,force=False):
    """Download Armadillo"""
    dep_download_dir=config.get('Main','dependency_download_dir')
    tools.download(arma_fname,arma_url,dep_download_dir,force)

def prepare(root,config):
    prefix=config.get('Main','prefix')
    dep_build_dir = config.get('Main','dependency_build_dir')
    dep_download_dir=config.get('Main','dependency_download_dir')
    arma_include_dir=prefix+"/bempp/include"

    print "Extracting Armadillo"
    extract_dir = dep_build_dir+"/"+arma_extract_dir
    tools.checkDeleteDirectory(extract_dir)
    try:
        tools.extract_file(dep_download_dir+"/"+arma_fname,dep_build_dir)
    except IOError:
        # Possibly a corrupted/truncated file. Try to download once again
        download(root,config,force=True)
        tools.extract_file(dep_download_dir+"/"+arma_fname,dep_build_dir)
    subprocess.check_call("cp -R "+extract_dir+"/include/* "+
                          prefix+"/bempp/include/",shell=True)
    print "Patching Armadillo"
    patch=py_patch.fromfile(root+"/installer/patches/armadillo_config.patch")
    cwd=os.getcwd()
    os.chdir(prefix+"/bempp/include/armadillo_bits")
    patch.apply()
    os.chdir(cwd)

    tools.setDefaultConfigOption(config,"Armadillo","include_dir",
                                 prefix+"/bempp/include",overwrite=True)

def configure(root,config):
    pass

def build(root,config):
    pass

def install(root,config):
    pass









