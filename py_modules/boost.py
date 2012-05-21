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


boost_fname='boost_1_49.tar.gz'
boost_url='http://www.bempp.org/files/boost_1_49.tar.gz'
boost_extract_dir='boost-cmake'
boost_dir='boost'
boost_version="1.49.0"

def configureBoost(root,config):
    """Download Boost if required"""

    boost_full_dir=root+"/contrib/"+boost_dir
    boost_download_name=root+"/contrib/files/"+boost_fname

    prefix=config.get('Main','prefix')
    boost_include_dir=prefix+"/bempp/contrib/boost/include/boost-"+boost_version
    
    if sys.platform.startswith('darwin'):
        unit_test_lib_name="libboost_unit_test_framework-mt.dylib"
    elif sys.platform.startswith('linux'):
        unit_test_lib_name="libboost_unit_test_framework-mt.so"
    else:
        raise Exception("Platform not supported")

    boost_unit_test_lib=prefix+"/bempp/contrib/boost/lib/boost-"+boost_version+"/"+unit_test_lib_name

    download_boost=True
    if config.has_option('Boost','download_boost'): download_boost=to_bool(config.get('Boost','download_boost'))
    
    if download_boost and not os.path.isdir(boost_full_dir):
        # Download boost

        if not os.path.isfile(boost_download_name):
            print "Downloading Boost ..."
            urllib.urlretrieve(boost_url,boost_download_name)

        print "Extracting Boost"
        extract_file(root+"/contrib/files/"+boost_fname,root+"/contrib/")
        os.rename(root+"/contrib/"+boost_extract_dir,boost_full_dir)
        shutil.copy(root+"/contrib/build_scripts/posix/boost_build.sh",boost_full_dir+"/boost_build.sh")

    if download_boost:
        if not config.has_section("Boost"): config.add_section("Boost")
        config.set("Boost",'unit_test_lib',boost_unit_test_lib)
        config.set("Boost",'include_dir',boost_include_dir)
        
    else:
        if not (config.has_option('Boost','unit_test_lib') and config.has_option('Boost','include_dir')):
            raise Exception("You need to specify 'unit_test_lib' and 'include_dir' under the 'Boost' header in the configuration file")

def buildBoost(root,config):

    download_boost=True
    if config.has_option('Boost','download_boost'): download_boost=to_bool(config.get('Boost','download_boost'))

    boost_full_dir=root+"/contrib/"+boost_dir

    prefix=config.get('Main','prefix')

    if download_boost and not os.path.isdir(prefix+"/bempp/contrib/boost"):
            print "Build Boost"
            cwd=os.getcwd()
            os.chdir(boost_full_dir)
            subprocess.call("sh ./boost_build.sh",shell=True)
            os.chdir(cwd)
    

        
    
        
        
        
            

        
