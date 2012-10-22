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

def download(root,config,force=False):
    dep_download_dir=config.get('Main','dependency_download_dir')
    swig_download_name=dep_download_dir+"/"+swig_fname
    tools.download(swig_fname,swig_url,dep_download_dir,force)

def prepare(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    dep_download_dir=config.get('Main','dependency_download_dir')
    prefix=config.get('Main','prefix')

    swig_full_dir=dep_build_dir+"/"+swig_dir
    swig_download_name=dep_download_dir+"/"+swig_fname

    swig_executable=prefix+"/bempp/bin/swig"

    tools.checkDeleteDirectory(swig_full_dir)

    print "Extracting Swig"
    try:
        tools.extract_file(dep_download_dir+"/"+swig_fname,dep_build_dir)
    except IOError:
        # Possibly a corrupted/truncated file. Try to download once again
        download(root,config,force=True)
        tools.extract_file(dep_download_dir+"/"+swig_fname,dep_build_dir)
    os.rename(dep_build_dir+"/"+swig_extract_dir,swig_full_dir)

    tools.setDefaultConfigOption(config,"Swig","exe",swig_executable,
                                 overwrite=True)
    tools.setCompilerOptions(config,'Swig')

def configure(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    swig_full_dir=dep_build_dir+"/"+swig_dir
    prefix=config.get('Main','prefix')
    swig_prefix=prefix+"/bempp"
    print "Configuring Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    tools.checkDeleteDirectory(swig_full_dir+"/build")
    subprocess.check_call(["./configure", "--prefix="+swig_prefix,
                           "--without-pcre", "--without-tcl",
                           "--without-perl5", "--without-octave",
                           "--without-java", "--without-gcj",
                           "--without-android", "--without-guile",
                           "--without-mzscheme", "--without-ruby",
                           "--without-ruby", "--without-php",
                           "--without-ocaml", "--without-pike",
                           "--without-chicken", "--without-csharp",
                           "--without-lua", "--without-allegrocl",
                           "--without-clisp", "--without-r",
                           "--without-go", "--without-d"])
    os.chdir(cwd)

def build(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    swig_full_dir=dep_build_dir+"/"+swig_dir
    njobs = tools.to_int(config.get('Main','build_jobs'))
    print "Building Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    subprocess.check_call(["make", "-j"+str(njobs)])
    os.chdir(cwd)

def install(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    swig_full_dir=dep_build_dir+"/"+swig_dir
    print "Installing Swig"
    cwd=os.getcwd()
    os.chdir(swig_full_dir)
    subprocess.check_call(["make", "install"])
    os.chdir(cwd)
