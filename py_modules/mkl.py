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
    lib_dir = prefix+"/bempp/lib"

    enable_mkl = tools.to_bool(tools.setDefaultConfigOption(config,'MKL',
                                                            'enable_mkl','no'))
    if enable_mkl:
        if config.has_option('BLAS', 'lib') and config.get('BLAS', 'lib'):
            print ("Warning: contents of the option 'lib' in section 'BLAS' "
                   "will be ignored, since enable_mkl is set")
        if config.has_option('LAPACK', 'lib') and config.get('LAPACK', 'lib'):
            print ("Warning: contents of the option 'lib' in section 'LAPACK' "
                   "will be ignored, since enable_mkl is set")
        
        # if mkl_source is not set, we'll get a sensible exception message
        mkl_source = config.get('MKL','mkl_source').lower()
        if mkl_source == 'installed':
            mkl_rt_lib = config.get('MKL','mkl_rt_lib')
            if not mkl_rt_lib:
                raise Exceptions("Option 'mkl_rt_lib' in section 'MKL' "
                                 "must not be empty")
            if not os.path.isfile(mkl_rt_lib):
                raise Exception("The path from option 'mkl_rt_lib' "
                                "in section 'MKL' is invalid")
            blas_lib = lapack_lib = mkl_rt_lib
        else:
	    mkl_link_files = ['libmkl_rt'] # MKL files that BEMPP links against.
            if mkl_source == 'redistributable':
                mkl_tarball=config.get('MKL','mkl_tarball')
                print 'Extracting MKL redistributables'
                tools.extract_file(mkl_tarball,prefix+"/bempp/lib/")
            elif mkl_source == 'enthought':
                print 'Creating symbolic links to Enthought MKL libraries'
                mkl_dir = sys.prefix+"/lib"
                if sys.platform.startswith('darwin'):
                    extension = ".dylib"
                elif sys.platform.startswith('linux'):
                    extension = ".so"
                else:
                    raise Exception("Unsupported platform: '"+sys.platform+"'")
                subprocess.check_call("ln -s "+mkl_dir+"/libmkl* "+lib_dir+"/",shell=True)
            else:
                raise Exception("Option 'mkl_source' in section 'MKL' must be "
                                "either 'installed', 'redistributable' or "
                                "'enthought'")
            blas_lib = ""
            if sys.platform.startswith('darwin'):
                for f in mkl_link_files: blas_lib = blas_lib+";"+lib_dir+"/"+f+".dylib"
            elif sys.platform.startswith('linux'):
                for f in mkl_link_files: blas_lib = blas_lib+";"+lib_dir+"/"+f+".so"
            else:
                raise Exception("Platform not supported")
            blas_lib = blas_lib[1:] # remove leading semicolon
            lapack_lib = blas_lib
            
        tools.setDefaultConfigOption(config,'BLAS','lib',blas_lib,overwrite=True)
        tools.setDefaultConfigOption(config,'LAPACK','lib',lapack_lib,overwrite=True)        
        
def configure(root,config):
    pass
             
def build(root,config):
    pass

def install(root,config):
    pass

        
    
        
        
        
            

        
