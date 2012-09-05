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

    debug=tools.setDefaultConfigOption(config,'Bempp','enable_debug','false')
    if tools.to_bool(debug):
        config.set('Bempp','build_type','Debug')
    else:
        config.set('Bempp','build_type','Release')

    tools.setDefaultConfigOption(config,'Bempp','build','true')

    tools.setDefaultConfigOption(config,'Bempp','build_dir',root+'/build')
    build_dir = config.get('Bempp','build_dir')
    # Replace ~ with /home/username
    build_dir = os.path.expanduser(build_dir)
    config.set('Bempp','build_dir',build_dir)

    tools.setCompilerOptions(config,'Bempp')

def configure(root,config):
    """Prepare the build of Bempp """

    tools.checkDeleteDirectory(config.get('Bempp','build_dir'))
    subprocess.check_call("sh ./.build.sh",shell=True)

def build(root,config):

    build_dir = config.get('Bempp','build_dir')
    do_build = tools.to_bool(config.get('Bempp','build'))
    njobs = tools.to_int(config.get('Main','build_jobs'))
    if do_build :
        cwd = os.getcwd()
        os.chdir(build_dir)
        subprocess.check_call("make -j"+str(njobs),shell=True)
        os.chdir(cwd)

def install(root,config):
    
    build_dir = config.get('Bempp','build_dir')
    do_build = tools.to_bool(config.get('Bempp','build'))
    if do_build :
        cwd = os.getcwd()
        os.chdir(build_dir)
        subprocess.check_call("make install",shell=True)
        os.chdir(cwd)
        
        
        
        
            

        
