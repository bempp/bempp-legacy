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
from py_modules.tools import extract_file, to_bool

import struct

tbb_fname_mac='tbb40_297oss_mac.tgz'
tbb_fname_linux='tbb40_297oss_lin.tgz'
tbb_url_mac='http://threadingbuildingblocks.org/uploads/78/181/4.0%20update%203/tbb40_297oss_mac.tgz'
tbb_url_linux='http://threadingbuildingblocks.org/uploads/78/181/4.0%20update%203/tbb40_297oss_lin.tgz'
tbb_extract_dir='tbb40_297oss'
tbb_dir='tbb'

def configureTbb(root,config):
    """Download Tbb if required"""

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

    prefix=config.get('Main','prefix')
    tbb_include_dir=prefix+"/bempp/contrib/tbb/include"
    
    if sys.platform.startswith('darwin'):
        tbb_lib_name="lib/libtbb.dylib"
        tbb_lib_name_debug="lib/libtbb_debug.dylib"
    elif sys.platform.startswith('linux'):
        if config.get('Main','architecture')=='intel64':
            # 64 bit
            tbb_lib_name="lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so"
            tbb_lib_name_debug="lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so"
        else:
            # presumably 32 bit
            tbb_lib_name="lib/ia32/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so"
            tbb_lib_name_debug="lib/ia32/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so"
    else:
        raise Exception("Platform not supported")


    download_tbb=True
    if config.has_option('Tbb','download_tbb'): download_tbb=to_bool(config.get('Tbb','download_tbb'))
    
    if download_tbb and not os.path.isdir(prefix+"/bempp/contrib/tbb"):
        # Download Tbb
        if not os.path.isfile(prefix+"/bempp/contrib/tbb"):
            print "Downloading Tbb ..."
            urllib.urlretrieve(tbb_url,root+"/contrib/files/"+tbb_fname)

        print "Extracting Tbb"
        extract_file(root+"/contrib/files/"+tbb_fname,prefix+"/bempp/contrib/")
        os.rename(prefix+"/bempp/contrib/"+tbb_extract_dir,prefix+"/bempp/contrib/"+tbb_dir)

    if download_tbb:
        if not config.has_section("Tbb"): config.add_section("Tbb")
        config.set("Tbb",'lib',prefix+"/bempp/contrib/tbb/"+tbb_lib_name)
        config.set("Tbb","lib_debug",prefix+"/bempp/contrib/tbb/"+tbb_lib_name_debug)
        config.set("Tbb",'include_dir',tbb_include_dir)        
    else:
        if not (config.has_option('Tbb','lib') and config.has_option('Tbb','lib_debug') and config.has_option('Tbb','include_dir')):
            raise Exception("You need to specify 'lib', 'lib_debug' and 'include_dir' under the 'Tbb' header in the configuration file")

     
def buildTbb(root,config):
    pass

    

        
    
        
        
        
            

        
