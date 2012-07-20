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


dune_fnames=['dune-common-2.1.0.tar.gz','dune-grid-2.1.0.tar.gz','dune-localfunctions-2.1.0.tar.gz']
dune_extract_names=['dune-common-2.1.0','dune-grid-2.1.0','dune-localfunctions-2.1.0']
dune_urls=['http://www.dune-project.org/download/2.1/dune-common-2.1.0.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-grid-2.1.0.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-localfunctions-2.1.0.tar.gz']
dune_names=['dune-common','dune-grid','dune-localfunctions']

def download(root,config):

    for i in range(3):
        tools.download(dune_fnames[i],dune_urls[i],root+"/contrib/files")

def prepare(root,config):

    prefix=config.get('Main','prefix')
    dune_dir=prefix+"/bempp/contrib/dune"

    # Download files
    for i in range(3):
        if not os.path.isdir(root+"/contrib/dune/"+dune_names[i]):
            print "Extract "+dune_names[i]
            tools.extract_file(root+"/contrib/files/"+dune_fnames[i],root+"/contrib/dune/")
            os.rename(root+"/contrib/dune/"+dune_extract_names[i],root+"/contrib/dune/"+dune_names[i])
            if i==1:
                print "Apply patch for "+dune_names[i]
                patch=py_patch.fromfile(root+"/contrib/patch/dune-grid.patch")
                cwd=os.getcwd()
                os.chdir(root+"/contrib/dune/dune-grid/dune/grid/utility")
                patch.apply()
                os.chdir(cwd)

def configure(root,config):
    prefix=config.get('Main','prefix')
    dune_dir=root+"/contrib/dune"
    dune_install_dir=prefix+"/bempp"
    cxx=config.get('Main','cxx')
    cc=config.get('Main','cc')
    print "Build Dune"
    cwd=os.getcwd()
    os.chdir(root+"/contrib/dune")
    f=open('dune_opts.ops','w')
    cflags = config.get('Main','cflags')
    cxxflags = config.get('Main','cxxflags')
    f.write("CONFIGURE_FLAGS=\" CXX="+cxx+" CC="+cc+" CXXFLAGS=\""+cxxflags+"\" CFLAGS=\""+cflags+"\" --enable-shared=yes --disable-documentation --enable-static=no --prefix="+dune_install_dir+"\"")
    f.close()
    subprocess.check_call("./dune-common/bin/dunecontrol --opts=./dune_opts.ops all",shell=True)
    os.chdir(cwd)

def build(root,config):
    cwd = os.getcwd()
    os.chdir(root+"/contrib/dune")
    subprocess.check_call("./dune-common/bin/dunecontrol make",shell=True)
    
def install(root,config):
    cwd = os.getcwd()
    ow.chdir(root+"/contrib/dune")
    subprocess.check_call("./dune-common/bin/dunecontrol make install",shell=True)
        
    
        
        
        
            

        
