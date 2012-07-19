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

import struct

mkl_fname='mkl.tar.gz'


def configureMkl(root,config):

    prefix=config.get('Main','prefix')
    lib_dir=prefix+"/bempp/lib"
    blas_files = ['libmkl_intel_lp64','libmkl_sequential','libmkl_core']
    lapack_files = ['libmkl_intel_lp64']

    if config.has_option('Main','enable_mkl'): enable_mkl=to_bool(config.get('Main','enable_mkl'))

    if enable_mkl:
        print "Extracting Mkl"
        extract_file(root+"/contrib/files/"+mkl_fname,prefix+"/bempp/lib/")

        blas_lib = ""
        lapack_lib = ""
        if sys.platform.startswith('darwin'):
            for f in blas_files: blas_lib = blas_lib+";"+lib_dir+"/"+f+".dylib"
            for f in lapack_files: lapack_lib = lapack_lib+";"+lib_dir+"/"+f+".dylib"                
        elif sys.platform.startswith('linux'):
            if config.get('Main','architecture')=='intel64':
                for f in blas_files: blas_lib = blas_lib+";"+lib_dir+"/intel64/"+f+".so"
                for f in lapack_files: lapack_lib = lapack_lib+";"+lib_dir+"/intel64/"+f+".so"                    
            elif config.get('Main','architecture')=='i386':
                for f in blas_files: blas_lib = blas_lib+";"+lib_dir+"/ia32/"+f+".so"            
                for f in lapack_files: lapack_lib = lapack_lib+";"+lib_dir+"/ia32/"+f+".so"                    
            else:
                raise Exception("Platform not supported")
        if not config.has_section('BLAS'): config.add_section('BLAS')
        if not config.has_section('LAPACK'): config.add_section('LAPACK')
        config.set('BLAS','lib',blas_lib)
        config.set('LAPACK','lib',lapack_lib)

        
     
def buildMkl(root,config):
    pass

    

        
    
        
        
        
            

        
