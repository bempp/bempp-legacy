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

def download(root,config):
    pass

def prepare(root,config):

    prefix = config.get('Main','prefix')
    lib_dir=prefix+"/bempp/lib"
    blas_files = ['libmkl_intel_lp64','libmkl_sequential','libmkl_core']
    lapack_files = ['libmkl_intel_lp64']

    enable_mkl = tools.to_bool(tools.setDefaultConfigOption(config,'Mkl','enable_mkl','no'))
    enthought_mkl = tools.to_bool(tools.setDefaultConfigOption(config,'Mkl','enable_enthought_mkl','no'))

    
    if enthought_mkl:
        print "Creating symbolic links to Enthought MKL Libraries"
        cwd = os.getcwd()
        mkl_dir = sys.prefix+"/lib"
        cmd_str = "ln -s "+mkl_dir+"/libmkl* ."
        os.chdir(lib_dir)
        subprocess.check_call(cmd_str,shell=True)
        os.chdir(cwd)
    elif enable_mkl:
        print "Extracting Mkl"
        if not config.has_option('Mkl','file_name'):
            raise Exception('Need to give full path of tar.gz archived file with Mkl redistributables')
        mkl_fname=config.get('Mkl','file_name')
        tools.extract_file(mkl_fname,prefix+"/bempp/lib/")

    if enable_mkl or enthought_mkl:
        blas_lib = ""
        lapack_lib = ""
        if sys.platform.startswith('darwin'):
            for f in blas_files: blas_lib = blas_lib+";"+lib_dir+"/"+f+".dylib"
            for f in lapack_files: lapack_lib = lapack_lib+";"+lib_dir+"/"+f+".dylib"                
        elif sys.platform.startswith('linux'):
            for f in blas_files: blas_lib = blas_lib+";"+lib_dir+"/"+f+".so"
            for f in lapack_files: lapack_lib = lapack_lib+";"+lib_dir+"/"+f+".so"                    
        else:
            raise Exception("Platform not supported")
        tools.setDefaultConfigOption(config,'BLAS','lib',blas_lib[1:],overwrite=True)
        tools.setDefaultConfigOption(config,'LAPACK','lib',lapack_lib[1:],overwrite=True)
        

def configure(root,config):
    pass
             
def build(root,config):
    pass

def install(root,config):
    pass

        
    
        
        
        
            

        
