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
from py_modules.tools import extract_file, to_bool, setDefaultConfigOption
from py_modules import python_patch as py_patch





def configureBempp(root,config):
    """Prepare the build of Bempp """

    debug=setDefaultConfigOption(config,'Bempp','enable_debug','false')
    if to_bool(debug):
        config.set('Bempp','build_type','Debug')
    else:
        config.set('Bempp','build_type','Release')

    setDefaultConfigOption(config,'Bempp','build','true')
    setDefaultConfigOption(config,'Bempp','build_dir',root+'/build')
    

    
def buildBempp(root,config):

    subprocess.call("sh .build.sh",shell=True)
    if to_bool(config.get('Bempp','build','true')):
        prefix=config.get('Main','prefix')
        build_dir=config.get('Bempp','build_dir')
        if not os.path.isdir(prefix+"/bempp/lib"):
            print "Build Bempp"
            subprocess.call("cd "+build_dir+"; make install",shell=True)

        
    
        
        
        
            

        
