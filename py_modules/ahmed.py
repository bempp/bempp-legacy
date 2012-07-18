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

ahmed_fname="AHMED-1.0.tar.gz"

def configureAhmed(root,config):

    if not config.has_section('AHMED'):
        config.add_section('AHMED')
        config.set('AHMED','enable_ahmed','OFF')
    else:
        if not config.has_option('AHMED','enable_ahmed'): 
            config.set('AHMED','enable_ahmed','ON')
    if to_bool(config.get('AHMED','enable_ahmed')):
        config.set('AHMED','enable_ahmed','ON') # Ensure that the option has the right format
        prefix=config.get('Main','prefix')
        ahmed_full_dir=root+"/contrib/ahmed"
        if os.path.isdir(ahmed_full_dir): shutil.rmtree(ahmed_full_dir)
        if sys.platform.startswith('darwin'):
            config.set('AHMED','lib',prefix+"/bempp/contrib/ahmed/lib/libAHMED.dylib")
        elif sys.platform.startswith('linux'):
            config.set('AHMED','lib',prefix+"/bempp/contrib/ahmed/lib/libAHMED.so")
        else:
            raise Exception("Platform not supported")
        config.set('AHMED','include_dir',prefix+"/bempp/contrib/ahmed/include/AHMED")
        print "Extracting AHMED"
        extract_file(root+"/contrib/files/"+ahmed_fname,root+"/contrib/")
        os.rename(root+"/contrib/AHMED_1.0",root+"/contrib/ahmed")
        shutil.copy(root+"/contrib/build_scripts/posix/ahmed_build.sh",trilinos_full_dir+"/ahmed_build.sh")

                 

def buildAhmed(root,config):

    if to_bool(config.get('AHMED','enable_ahmed')):
        print "Build AHMED"
        cwd=os.getcwd()
        os.chdir(root+"/contrib/ahmed")
        subprocess.check_call("sh ./ahmed_build.sh",shell=True)
        os.chdir(cwd)


        
    
        
        
        
            

        
