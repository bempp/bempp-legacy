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


dune_fnames=['dune-common-2.1.0.tar.gz','dune-grid-2.1.0.tar.gz','dune-localfunctions-2.1.0.tar.gz']
dune_extract_names=['dune-common-2.1.0','dune-grid-2.1.0','dune-localfunctions-2.1.0']
dune_urls=['http://www.dune-project.org/download/2.1/dune-common-2.1.0.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-grid-2.1.0.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-localfunctions-2.1.0.tar.gz']
dune_names=['dune-common','dune-grid','dune-localfunctions']



def checkAndBuildDune(root,config):
    """Download and build Dune if required"""


    prefix=config.get('Main','prefix')
    dune_dir=prefix+"/bempp/contrib/dune"
    
    

    download_dune=True
    if not os.path.isdir(prefix+"/bempp/contrib/dune"):
        os.mkdir(dune_dir)
        # Download files
        for i in range(3):
            if not os.path.isfile(root+"/contrib/files/"+dune_fnames[i]):
                print "Download "+dune_names[i]
                urllib.urlretrieve(dune_urls[i],root+"/contrib/files/"+dune_fnames[i])
        # Extract and patch
        for i in range(3):
            print "Extract "+dune_names[i]
            extract_file(root+"/contrib/files/"+dune_fnames[i],prefix+"/bempp/contrib/dune/")
            os.rename(prefix+"/bempp/contrib/dune/"+dune_extract_names[i],prefix+"/bempp/contrib/dune/"+dune_names[i])
        shutil.copytree(root+"/contrib/dune/dune-foamgrid",prefix+"/bempp/contrib/dune/dune-foamgrid")
        print "Apply patch"
        patch=py_patch.fromfile(root+"/contrib/patch/dune-grid.patch")
        cwd=os.getcwd()
        os.chdir(prefix+"/bempp/contrib/dune/dune-grid/dune/grid/utility")
        patch.apply()
        os.chdir(cwd)
        print "Build Dune"
        os.chdir(prefix+"/bempp/contrib/dune")
        subprocess.call("./dune-common/bin/dunecontrol all",shell=True)

    

        
    
        
        
        
            

        
