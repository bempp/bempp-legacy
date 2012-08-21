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

import struct

tbb_fname_mac='tbb40_297oss_mac.tgz'
tbb_fname_linux='tbb40_297oss_lin.tgz'
tbb_url_mac='http://threadingbuildingblocks.org/uploads/78/181/4.0%20update%203/tbb40_297oss_mac.tgz'
tbb_url_linux='http://threadingbuildingblocks.org/uploads/78/181/4.0%20update%203/tbb40_297oss_lin.tgz'
tbb_extract_dir='tbb40_297oss'
tbb_dir='tbb'
tbb_fname_short='tbb.tgz'

def download(root,config):
    tbb_full_dir=root+"/contrib/"+tbb_dir
    if sys.platform.startswith('darwin'):
        tbb_download_name=root+"/contrib/files/"+tbb_fname_mac
        tbb_url=tbb_url_mac
        tbb_fname=tbb_fname_mac
    elif sys.platform.startswith('linux'):
        tbb_download_name=root+"/contrib/files/"+tbb_fname_linux
        tbb_url=tbb_url_linux
        tbb_fname=tbb_fname_linux
    else:
        raise Exception("Platform not supported")
    tools.download(tbb_fname_short,tbb_url,root+"/contrib/files")

def prepare(root,config):

    prefix=config.get('Main','prefix')

    print "Extracting Tbb"

    tools.checkDeleteDirectory(root+"/contrib/tbb")
    tools.extract_file(root+"/contrib/files/"+tbb_fname_short,root+"/contrib/")
    os.rename(root+"/contrib/"+tbb_extract_dir,root+"/contrib/tbb")
    subprocess.check_call("cp -R "+root+"/contrib/tbb/include/* "+prefix+"/bempp/include/",shell=True)


    if sys.platform.startswith('darwin'):
        libdir_orig = root+"/contrib/tbb/lib"
        tbb_lib_name="libtbb.dylib"
        tbb_lib_name_debug="libtbb_debug.dylib"
    elif sys.platform.startswith('linux'):
        tbb_lib_name = "libtbb.so"
        tbb_lib_name_debug = "libtbb_debug.so"
        arch = config.get('Main','architecture')
        if arch in ('intel64','ia32','ia64'):
            libdir_orig = root+"/contrib/tbb/lib/"+arch+"/cc4.1.0_libc2.4_kernel2.6.16.21"
        else:
            raise Exception("Unrecognized architecture: '"+arch+"'")
    else:
        raise Exception("Platform not supported")

    subprocess.check_call("cp -R "+libdir_orig+"/* "+prefix+"/bempp/lib/",shell=True)

    tools.setDefaultConfigOption(config,"Tbb",'lib',prefix+"/bempp/lib/"+tbb_lib_name,overwrite=True)
    tools.setDefaultConfigOption(config,"Tbb","lib_debug",prefix+"/bempp/lib/"+tbb_lib_name_debug,overwrite=True)
    tools.setDefaultConfigOption(config,"Tbb",'include_dir',prefix+"/bempp/include",overwrite=True) 
    

def configure(root,config):
    pass

def build(root,config):
    pass

def install(root,config):
    pass
    

        
    
        
        
        
            

        
