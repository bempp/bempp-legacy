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

import glob,os,urllib,shutil,subprocess,sys
from py_modules import tools 

def find_file_in_dirs(fname, dirs):
    for d in dirs:
        path = os.path.join(d, fname)
        if os.path.isfile(path):
            return path
    return None

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

        mkl_fname_cores = ["mkl_intel_thread",
                           "mkl_core","iomp5"]
        arch = config.get('Main','architecture')
        if arch == 'ia32':
            mkl_fname_cores.append("mkl_intel")
        else:
            mkl_fname_cores.append("mkl_intel_lp64")
        
        # If mkl_source is not set, we'll get a sensible exception message
        mkl_source = config.get('MKL','mkl_source').lower()

        # Set the variable 'mkl_dirs' to a list of directories
        # where mkl libraries are to be found
        if mkl_source == 'installed':
            mkl_dirs = config.get('MKL','mkl_dirs')
            mkl_dirs = mkl_dirs.split(";")
            mkl_dirs = [s.strip() for s in mkl_dirs]
            if not mkl_dirs:
                raise Exception("Option 'mkl_dirs' in section 'MKL' "
                                "must not be empty")
            for d in mkl_dirs:
                if not os.path.isdir(d):
                    raise Exception("Option 'mkl_dirs' in section 'MKL' "
                                    "is invalid: '"+d+"' is not a directory")
        elif mkl_source == 'redistributable':
            mkl_tarball = config.get('MKL','mkl_tarball')
            print 'Extracting MKL redistributables'
            tools.extract_file(mkl_tarball,lib_dir)
            mkl_dirs = [lib_dir]
        elif mkl_source == 'enthought':
            d = os.path.join(sys.prefix,"lib")
            if not os.path.isdir(d):
                raise Exception("Unexpected error: '"+d+"' is not a directory")
            mkl_dirs = [d]
        else:
            raise Exception("Option 'mkl_source' in section 'MKL' must be "
                            "either 'installed', 'redistributable' or "
                            "'enthought'")

        # Find MKL library files and create symbolic links to them
        # in the lib_dir directory
        if sys.platform.startswith('darwin'):
            extension = ".dylib"
        elif sys.platform.startswith('linux'):
            extension = ".so"
        else:
            raise Exception("Unsupported platform: '"+sys.platform+"'")
        if mkl_source != 'redistributable':
            for f in mkl_fname_cores:
                fname = "lib"+f+extension
                path = find_file_in_dirs(fname, mkl_dirs)
                if not path:
                    raise Exception("MKL library '"+fname+"' not found in any"
                                    "of these directories: '"+
                                    "', '".join(mkl_dirs)+"'")
                os.symlink(path, os.path.join(lib_dir,fname))

            # Symlink to all other libmkl* files found in mkl_dir directories
            for d in mkl_dirs:
                for path in glob.glob(os.path.join(d,"libmkl")+"*"+extension):
                    fname = os.path.basename(path)
                    new_path = os.path.join(lib_dir,fname)
                    if not os.path.exists(new_path):
                        os.symlink(path, new_path)

        blas_lib = ""
        for f in mkl_fname_cores:
            blas_lib = blas_lib + os.path.join(lib_dir,"lib"+f+extension)+";"
        blas_lib = blas_lib+";-pthread"
        lapack_lib = blas_lib

        tools.setDefaultConfigOption(config,'BLAS','lib',blas_lib,overwrite=True)
        tools.setDefaultConfigOption(config,'LAPACK','lib',lapack_lib,overwrite=True)
    else: # enable_mkl is false
        # Create symbolic links to BLAS and LAPACK libraries in
        # lib_dir and replace original locations in the linking line
        # with these symbolic links
        for section in ('BLAS','LAPACK'):
            paths = config.get(section,'lib')
            paths = paths.split(";")
            paths = [s.strip() for s in paths]
            new_paths = []
            for path in paths:
                if os.path.isfile(path):
                    fname = os.path.basename(path)
                    new_path = os.path.join(lib_dir,fname)
                    if not os.path.exists(new_path):
                        os.symlink(path, new_path)
                    new_paths.append(new_path)
                else:
                    new_paths.append(path)
            new_setting = ";".join(new_paths)
            tools.setDefaultConfigOption(config,section,'lib',
                                         new_setting,overwrite=True)
        
def configure(root,config):
    pass
             
def build(root,config):
    pass

def install(root,config):
    pass

        
    
        
        
        
            

        
