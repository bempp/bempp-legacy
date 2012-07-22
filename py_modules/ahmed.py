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

def download(root,config):
    pass

def prepare(root,config):

    enable_ahmed = tools.to_bool(tools.setDefaultConfigOption(config,'AHMED','enable_ahmed','false'))
    if enable_ahmed:
        if not config.has_option('AHMED','file_name'):
            raise Exception('Need to give full path of tar.gz archived file with AHMED 1.0 release')
        ahmed_fname=config.get('AHMED','file_name')
        config.set('AHMED','with_ahmed','ON') 
        prefix=config.get('Main','prefix')
        arch = config.get('Main','architecture')
        if arch == 'i386':
            config.set('AHMED','enable64','OFF')
        else:
            config.set('AHMED','enable64','ON')
        ahmed_full_dir=root+"/contrib/ahmed"
        tools.checkDeleteDirectory(ahmed_full_dir)
        if sys.platform.startswith('darwin'):
            config.set('AHMED','lib',prefix+"/bempp/contrib/ahmed/lib/libAHMED.dylib")
        elif sys.platform.startswith('linux'):
            config.set('AHMED','lib',prefix+"/bempp/contrib/ahmed/lib/libAHMED.so")
        else:
            raise Exception("Platform not supported")
        config.set('AHMED','include_dir',prefix+"/bempp/include/AHMED")
        print "Extracting AHMED"
        tools.extract_file(ahmed_fname,root+"/contrib/")
        os.rename(root+"/contrib/AHMED_1.0",root+"/contrib/ahmed")
        shutil.copy(root+"/contrib/build_scripts/posix/ahmed_build.sh",ahmed_full_dir+"/ahmed_build.sh")
    else:
        config.set('AHMED','with_ahmed','OFF')

def configure(root,config):
    prefix = config.get('Main','prefix')
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        print "Configure AHMED"
        cwd=os.getcwd()
        os.chdir(root+"/contrib/ahmed")
        tools.checkDeleteDirectory(root+"/contrib/ahmed/build")
        subprocess.check_call("sh ./ahmed_build.sh",shell=True)
        os.chdir(cwd)
    

def build(root,config):
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        njobs = tools.to_int(config.get('Main','build_jobs',1))
        print "Build AHMED"
        cwd=os.getcwd()
        os.chdir(root+"/contrib/ahmed/build")
        subprocess.check_call("make -j"+str(njobs),shell=True)
        os.chdir(cwd)

def install(root,config):
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        prefix = config.get('Main','prefix')
        print "Install AHMED"
        cwd=os.getcwd()
        os.chdir(root+"/contrib/ahmed/build")
        subprocess.check_call("make install",shell=True)
        os.chdir(prefix+"/bempp/include/AHMED")
        g77 = tools.to_bool(config.get('AHMED','with_g77'))
        print "Patching AHMED for G77 calling BLAS convention"
        patch=py_patch.fromfile(root+"/contrib/patch/ahmed_blas.patch")
        patch.apply()
        os.chdir(cwd)
        
        
        
        
            

        
