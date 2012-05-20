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


trilinos_fname='trilinos-10.10.2-Source.tar.gz'
trilinos_extract_dir='trilinos-10.10.2-Source'
trilinos_dir='trilinos'
trilinos_url='http://trilinos.sandia.gov/download/files/trilinos-10.10.2-Source.tar.gz'



def checkAndBuildTrilinos(root,config):
    """Download and build Trilinos if required"""

    trilinos_full_dir=root+"/contrib/"+trilinos_dir
    trilinos_download_name=root+"/contrib/files/"+trilinos_fname

    prefix=config.get('Main','prefix')    
    download_trilinos=True
    if config.has_option('Trilinos','download_trilinos'): download_trilinos=to_bool(config.get('Trilinos','download_trilinos'))
    
    if download_trilinos and not os.path.isdir(trilinos_full_dir):
        # Download Trilinos

        if not os.path.isfile(trilinos_download_name):
            print "Downloading Trilinos ..."
            urllib.urlretrieve(trilinos_url,trilinos_download_name)

        print "Extracting Trilinos"
        extract_file(root+"/contrib/files/"+trilinos_fname,root+"/contrib/")
        os.rename(root+"/contrib/"+trilinos_extract_dir,root+"/contrib/"+trilinos_dir)
        shutil.copy(root+"/contrib/build_scripts/posix/trilinos_build.sh",trilinos_full_dir+"/trilinos_build.sh")
        print "Patching ..."
        patch=py_patch.fromfile(root+"/contrib/patch/Thyra_BelosLinearOpWithSolve_def.patch")
        cwd=os.getcwd()
        os.chdir(root+"/contrib/trilinos//packages/stratimikos/adapters/belos/src")
        patch.apply()
        os.chdir(cwd)
        

    if download_trilinos:
        if not os.path.isdir(prefix+"/bempp/contrib/trilinos"):
            print "Install Trilinos"
            cwd=os.getcwd()
            os.chdir(trilinos_full_dir)
            subprocess.call("sh ./trilinos_build.sh",shell=True)
        if not config.has_section("Trilinos"): config.add_section("Trilinos")
        config.set("Trilinos",'cmake_path',prefix+"/bempp/contrib/trilinos/lib/cmake/Trilinos/")
        
    else:
        if not config.has_option('Trilinos','cmake_path'):
            raise Exception("You need to specify 'cmake_path' under the 'Trilinos' header in the configuration file")

    

        
    
        
        
        
            

        
