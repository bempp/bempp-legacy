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


swig_fname='swig-2.0.8.tar.gz'
swig_url='http://prdownloads.sourceforge.net/swig/swig-2.0.8.tar.gz'
swig_extract_dir='swig-2.0.8'
swig_dir='swig'

def download(root,config):
    swig_download_name=root+"/contrib/files/"+swig_fname
    tools.download(swig_fname,swig_url,root+"/contrib/files")

def prepare(root,config):
    swig_full_dir=root+"/contrib/"+swig_dir
    swig_download_name=root+"/contrib/files/"+swig_fname

    prefix=config.get('Main','prefix')
    swig_executable=prefix+"/bempp/bin/swig"

    tools.checkDeleteDirectory(swig_full_dir)

    print "Extracting Swig"
    tools.extract_file(root+"/contrib/files/"+swig_fname,root+"/contrib/")
    os.rename(root+"/contrib/"+swig_extract_dir,swig_full_dir)

    tools.setDefaultConfigOption(config,"Swig","exe",swig_executable,
                                 overwrite=True)
    tools.setCompilerOptions(config,'Swig')

def configure(root,config):
    swig_full_dir=root+"/contrib/"+swig_dir
    prefix=config.get('Main','prefix')
    swig_prefix=prefix+"/bempp"
    print "Configuring Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    tools.checkDeleteDirectory(swig_full_dir+"/build")
    subprocess.check_call(["./configure", "--prefix="+swig_prefix, "--without-pcre"])
    os.chdir(cwd)

def build(root,config):
    swig_full_dir=root+"/contrib/"+swig_dir
    njobs = tools.to_int(config.get('Main','build_jobs'))
    print "Build Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    subprocess.check_call(["make", "-j"+str(njobs)])
    os.chdir(cwd)

def install(root,config):
    swig_full_dir=root+"/contrib/"+swig_dir
    print "Install Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    subprocess.check_call(["make", "install"])
    os.chdir(cwd)
