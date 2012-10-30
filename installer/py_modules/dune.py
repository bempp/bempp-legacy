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


dune_fnames=['dune-common-2.1.1.tar.gz',
             'dune-grid-2.1.1.tar.gz',
             'dune-localfunctions-2.1.1.tar.gz']
dune_extract_names=['dune-common-2.1.1',
                    'dune-grid-2.1.1',
                    'dune-localfunctions-2.1.1']
dune_urls=['http://www.dune-project.org/download/2.1/dune-common-2.1.1.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-grid-2.1.1.tar.gz',
	   'http://www.dune-project.org/download/2.1/dune-localfunctions-2.1.1.tar.gz']
dune_names=['dune-common','dune-grid','dune-localfunctions']

def download(root,config,force=False):
    dep_download_dir=config.get('Main','dependency_download_dir')
    for i in range(3):
        tools.download(dune_fnames[i],dune_urls[i],dep_download_dir,force)

def prepare(root,config):
    dep_build_dir=config.get('Main','dependency_build_dir')
    dep_download_dir=config.get('Main','dependency_download_dir')

    # Download files
    tools.checkCreateDir(dep_build_dir+"/dune")
    for i in range(3):
        tools.checkDeleteDirectory(dep_build_dir+"/dune/"+dune_names[i])
        print "Extracting "+dune_names[i]
        try:
            tools.extract_file(dep_download_dir+"/"+dune_fnames[i],
                               dep_build_dir+"/dune/")
        except IOError:
            # Possibly a corrupted/truncated file. Try to download once again
            download(root,config,force=True)
            tools.extract_file(dep_download_dir+"/"+dune_fnames[i],
                               dep_build_dir+"/dune/")
        os.rename(dep_build_dir+"/dune/"+dune_extract_names[i],
                  dep_build_dir+"/dune/"+dune_names[i])
        if i==1:
            print "Patching "+dune_names[i]
            patch=py_patch.fromfile(root+"/installer/patches/dune-grid.patch")
            cwd=os.getcwd()
            os.chdir(dep_build_dir+"/dune/dune-grid/dune/grid/utility")
            patch.apply()
            os.chdir(cwd)
    # Copy over dune-foamgrid
    tools.checkDeleteDirectory(dep_build_dir+"/dune/dune-foamgrid")
    shutil.copytree(root+"/contrib/dune/dune-foamgrid",
                    dep_build_dir+"/dune/dune-foamgrid")

    tools.setCompilerOptions(config,'Dune')

def configure(root,config):
    prefix=config.get('Main','prefix')
    dep_build_dir=config.get('Main','dependency_build_dir')
    dune_dir=dep_build_dir+"/dune"
    dune_install_dir=prefix+"/bempp"
    cxx=config.get('Dune','cxx')
    cc=config.get('Dune','cc')
    cflags = config.get('Dune','cflags')
    cxxflags = config.get('Dune','cxxflags')
    cwd=os.getcwd()

    njobs = tools.to_int(config.get('Main','build_jobs'))
    config_string_common = ("CC="+cc+" CXX="+cxx+" CFLAGS='"+cflags+
                            "' CXXFLAGS='"+cxxflags+
                            "' ./configure --disable-documentation "
                            "--enable-shared=yes --enable-static=no "
                            "--disable-gxx0xcheck "
                            "--prefix='"+dune_install_dir+"' "
                            "--libdir='"+dune_install_dir+"/lib'")
    config_string_grid = (config_string_common+
                          " --with-dune-common="+dune_dir+"/dune-common")
    config_string_localfunctions = (config_string_grid+
                                    " --with-dune-grid="+dune_dir+"/dune-grid")
    os.chdir(dune_dir+"/dune-common")
    subprocess.check_call(config_string_common,shell=True)
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

    os.chdir(dune_dir+"/dune-grid")
    subprocess.check_call(config_string_grid,shell=True)
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

    os.chdir(dune_dir+"/dune-localfunctions")
    subprocess.check_call(config_string_localfunctions,shell=True)
    subprocess.check_call("make -j"+str(njobs),shell=True)
    os.chdir(cwd)

def build(root,config):
    pass # Has already been built in the configure step


def install(root,config):
    prefix=config.get('Main','prefix')
    dep_build_dir=config.get('Main','dependency_build_dir')
    dune_dir=dep_build_dir+"/dune"
    cwd=os.getcwd()

    os.chdir(dune_dir+"/dune-common")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)

    os.chdir(dune_dir+"/dune-grid")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)

    os.chdir(dune_dir+"/dune-localfunctions")
    subprocess.check_call("make install",shell=True)
    os.chdir(cwd)

    os.chdir(dune_dir+"/dune-foamgrid")
    subprocess.check_call("find . -name '*.hh' | cpio -pdm "+
                          prefix+"/bempp/include",shell=True)
    os.chdir(cwd)







