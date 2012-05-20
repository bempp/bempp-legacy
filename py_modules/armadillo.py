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
from py_modules import python_patch as py_patch


arma_fname='armadillo-3.2.0.tar.gz'
arma_url='http://sourceforge.net/projects/arma/files/armadillo-3.2.0.tar.gz'
arma_extract_dir='armadillo-3.2.0'
arma_dir='armadillo'

def checkAndBuildArmadillo(root,config):
    """Download and build Armadillo if required"""


    prefix=config.get('Main','prefix')
    arma_include_dir=prefix+"/bempp/contrib/armadillo/include"
    

    download_arma=True
    if config.has_option('Armadillo','download_armadillo'): download_arma=to_bool(config.get('Armadillo','download_armadillo'))
    
    if download_arma and not os.path.isdir(prefix+"/bempp/contrib/armadillo"):
        # Download Armadillo
        if not os.path.isfile(root+"/contrib/files/"+arma_fname):
            print "Downloading Armadillo ..."
            urllib.urlretrieve(arma_url,root+"/contrib/files/"+arma_fname)

        print "Extracting Armadillo"
        extract_file(root+"/contrib/files/"+arma_fname,prefix+"/bempp/contrib/")
        os.rename(prefix+"/bempp/contrib/"+arma_extract_dir,prefix+"/bempp/contrib/"+arma_dir)
        print "Applying patches"
        patch=py_patch.fromfile(root+"/contrib/patch/armadillo_config.patch")
        cwd=os.getcwd()
        os.chdir(prefix+"/bempp/contrib/armadillo/include/armadillo_bits")
        patch.apply()
        os.chdir(cwd)
        

    if download_arma:
        if not config.has_section("Armadillo"): config.add_section("Armadillo")
        config.set("Armadillo",'include_dir',prefix+"/bempp/contrib/armadillo/include")
    else:
        if not config.has_option('Armadillo','include_dir'):
            raise Exception("You need to specify 'include_dir' under the 'Armadillo' header in the configuration file")

     

    

        
    
        
        
        
            

        
