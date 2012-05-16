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

# This script downloads necessary third party libraries, which are needed for BEM++


import sys
import os
import urllib
import tarfile
import zipfile
import shutil

# Files to download

file_list={'trilinos':('trilinos-10.10.2-Source.tar.gz','http://trilinos.sandia.gov/download/files/trilinos-10.10.2-Source.tar.gz'),
	   'armadillo':('armadillo-3.0.3.tar.gz','http://sourceforge.net/projects/arma/files/armadillo-3.0.3.tar.gz'),
           'dune-common':('dune-common-2.1.0.tar.gz','http://www.dune-project.org/download/2.1/dune-common-2.1.0.tar.gz'),
	   'dune-grid':('dune-grid-2.1.0.tar.gz','http://www.dune-project.org/download/2.1/dune-grid-2.1.0.tar.gz'),
	   'dune-localfunctions':('dune-localfunctions-2.1.0.tar.gz','http://www.dune-project.org/download/2.1/dune-localfunctions-2.1.0.tar.gz'),
           'boost':('boost_1_49.tar.gz','http://www.bempp.org/files/boost_1_49.tar.gz')}

rename_list={'trilinos-10.10.2-Source':'trilinos',
	     'armadillo-3.0.3':'armadillo',
	     'dune-common-2.1.0':'dune/dune-common',
	     'dune-grid-2.1.0':'dune/dune-grid',
	     'dune-localfunctions-2.1.0':'dune/dune-localfunctions'}

patch_list=[('dune-grid.patch','dune/dune-grid/dune/grid/utility'),
	    ('Thyra_BelosLinearOpWithSolve_def.patch','trilinos/packages/stratimikos/adapters/belos/src'),
	    ('armadillo_config.patch','armadillo/include/armadillo_bits')]	

###########################

# The following functions determine the directory, in which the file resides.
# Trick taken from http://stackoverflow.com/questions/2632199/how-do-i-get-the-path-of-the-current-executed-file-in-python

def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def module_path():
    if we_are_frozen():
        return os.path.split(os.path.abspath(sys.executable))[0]
    return os.path.split(os.path.abspath(__file__))[0]

###########################

def download_files():
    # Check if files exist and if not download them.
    fpath=module_path()+"/contrib/files/"
    for item in file_list:
	    fname=fpath+file_list[item][0]
	    if not os.path.isfile(fname): 
	    	    print "Downloading "+file_list[item][0]+" ..."
		    urllib.urlretrieve(file_list[item][1],fname)


##########################

# The following function extracts the tar gz files. It is taken from
# http://code.activestate.com/recipes/576714-extract-a-compressed-file/

def extract_file(path, to_directory='.'):
    if path.endswith('.zip'):
        opener, mode = zipfile.ZipFile, 'r'
    elif path.endswith('.tar.gz') or path.endswith('.tgz'):
        opener, mode = tarfile.open, 'r:gz'
    elif path.endswith('.tar.bz2') or path.endswith('.tbz'):
        opener, mode = tarfile.open, 'r:bz2'
    else: 
        raise ValueError, "Could not extract `%s` as no appropriate extractor is found" % path
    
    cwd = os.getcwd()
    os.chdir(to_directory)
    
    try:
        file = opener(path, mode)
        try: file.extractall()
        finally: file.close()
    finally:
        os.chdir(cwd)

##########################

def clean_up():
   # Clean up existing files in contrib directory
   fpath=module_path()+"/contrib/"
   for item in rename_list:
	   name=fpath+rename_list[item]
   	   if os.path.isdir(name):
	    	print rename_list[item]+" exists, deleting ..."
	    	shutil.rmtree(name)
	

def extract_all_files():
    # Extract all downloaded files
    fpath=module_path()+"/contrib/files/"
    extract_dir=module_path()+"/contrib/"
    for item in file_list:
	    print "Extracting "+file_list[item][0]
	    fname=fpath+file_list[item][0]
	    extract_file(fname,extract_dir)

def rename_directories():
    # Strip version numbers from download directory names
    fpath=module_path()+"/contrib/"
    for item in rename_list:
	    name=fpath+rename_list[item]
    	    os.rename(fpath+item,name)

def apply_patches():
    import python_patch as p
    fpath=module_path()+"/contrib/"
    cwd=os.getcwd()
    print "Apply patches ..."
    for item in patch_list:
	    patch=p.fromfile(fpath+"patch/"+item[0])
	    os.chdir(fpath+item[1])
	    patch.apply()
	    os.chdir(cwd)
	

def install_cmake_files():
    dune_cmake_files()

def check_platform_for_tbb():
    if sys.platform=="darwin":
        file_list['tbb']=('tbb40_20120408oss_mac.tgz','http://threadingbuildingblocks.org/uploads/77/185/4.0%20update%204/tbb40_20120408oss_mac.tgz')
        rename_list['tbb']=('tbb40_20120408oss','tbb')
    elif sys.platform=='linux2':
        file_list['tbb']=('tbb40_20120408oss_lin.tgz','http://threadingbuildingblocks.org/uploads/77/185/4.0%20update%204/tbb40_20120408oss_lin.tgz')
        rename_list['tbb']=('tbb40_20120408oss','tbb')
    else:
        raise Exception("Platform not supported")

def dune_cmake_files():
    src=module_path()+"/contrib/patch/cmake/"
    dst=module_path()+"/contrib"
    shutil.copy(src+"dune-common.cmake",dst+"/dune/dune-common/CMakeLists.txt")
    shutil.copy(src+"dune-grid.cmake",dst+"/dune/dune-grid/CMakeLists.txt")
    shutil.copy(src+"dune-foamgrid.cmake",dst+"/dune/dune-foamgrid/CMakeLists.txt")
    shutil.copy(src+"dune-localfunctions.cmake",dst+"/dune/dune-localfunctions/CMakeLists.txt")


###########################

if __name__ == "__main__":
        check_platform_for_tbb()
	download_files()
	clean_up()
	extract_all_files()
	rename_directories()
	apply_patches()
	install_cmake_files()
        print "Please enter the subdirectory contrib and run \n cmake . && make install \n to install the dependencies"


