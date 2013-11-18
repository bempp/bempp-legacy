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
        ahmed_fname=tools.normalizePath(config,ahmed_fname)
        config.set('AHMED','with_ahmed','ON')
        prefix=config.get('Main','prefix')
        dep_build_dir = config.get('Main','dependency_build_dir')
        arch = config.get('Main','architecture')
        if arch == 'ia32':
            config.set('AHMED','enable64','OFF')
        else:
            config.set('AHMED','enable64','ON')
        ahmed_full_dir=dep_build_dir+"/ahmed"
        tools.checkDeleteDirectory(ahmed_full_dir)
        if sys.platform.startswith('darwin'):
            config.set('AHMED','lib',prefix+"/bempp/lib/libAHMED.dylib")
        elif sys.platform.startswith('linux'):
            config.set('AHMED','lib',prefix+"/bempp/lib/libAHMED.so")
        else:
            raise Exception("Platform not supported")
        config.set('AHMED','include_dir',prefix+"/bempp/include/AHMED")
        if not os.path.isfile(ahmed_fname):
            raise Exception("File '"+ahmed_fname+"' does not exist")
        print "Extracting AHMED"
        tools.extract_file(ahmed_fname,dep_build_dir)
        os.rename(dep_build_dir+"/AHMED_1.0",ahmed_full_dir)
        shutil.copy(root+"/installer/build_scripts/posix/ahmed_build.sh",ahmed_full_dir+"/ahmed_build.sh")
        print "Patching AHMED"
        cwd=os.getcwd()
        os.chdir(ahmed_full_dir)
        for s in ("ahmed_cmake.patch",
                  "ahmed_addGeHGeH_single_precision.patch",
                  "ahmed_changelog_H.h.patch",
                  "ahmed_pass_clusters_to_aca_matgen_apprx.h.patch",
                  "ahmed_generic_aca_apprx.h.patch",
                  "ahmed_check_error_apprx.h.patch",
                  "ahmed_bbx_apprx.h.patch",
                  "ahmed_pass_clusters_to_aca_matgen_ACA.h.patch",
                  "ahmed_frobenius_norm_ACA.h.patch",
                  "ahmed_zero_pu_pv_ACA.h.patch",
                  "ahmed_retry_if_zero_and_orig_cross_ACA.h.patch",
                  "ahmed_changelog_ACA.h.patch",
                  "ahmed_permuted_indices_bemcluster.h.patch",
                  "ahmed_changelog_bemcluster.h.patch",
		  "ahmed_omp.patch",
                  "ahmed_basmod_h.patch"):
            py_patch.fromfile(root+"/installer/patches/"+s).apply()
        shutil.copy(root+"/installer/patches/ahmed_bbx_bbxbemcluster.h",
                    "./Include/bbxbemcluster.h")
        shutil.copy(root+"/installer/patches/ahmed_bbx_bbxbemblcluster.h",
                    "./Include/bbxbemblcluster.h")
        os.chdir(cwd)
        tools.setCompilerOptions(config,'AHMED')
    else:
        config.set('AHMED','with_ahmed','OFF')


def configure(root,config):
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        dep_build_dir = config.get('Main','dependency_build_dir')
        print "Configuring AHMED"
        cwd=os.getcwd()
        os.chdir(dep_build_dir+"/ahmed")
        tools.checkDeleteDirectory(dep_build_dir+"/ahmed/build")
        subprocess.check_call("sh ./ahmed_build.sh",shell=True)
        os.chdir(cwd)


def build(root,config):
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        dep_build_dir = config.get('Main','dependency_build_dir')
        njobs = tools.to_int(config.get('Main','build_jobs'))
        print "Building AHMED"
        cwd=os.getcwd()
        os.chdir(dep_build_dir+"/ahmed/build")
        subprocess.check_call("make -j"+str(njobs),shell=True)
        os.chdir(cwd)

def install(root,config):
    if tools.to_bool(config.get('AHMED','enable_ahmed')):
        dep_build_dir = config.get('Main','dependency_build_dir')
        prefix = config.get('Main','prefix')
        print "Installing AHMED"
        cwd=os.getcwd()
        os.chdir(dep_build_dir+"/ahmed/build")
        subprocess.check_call("make install",shell=True)
        g77 = tools.to_bool(config.get('AHMED','with_g77'))
        if g77:
            os.chdir(prefix+"/bempp/include/AHMED")
            print "Patching AHMED for G77 calling BLAS convention"
            patch=py_patch.fromfile(root+"/installer/patches/ahmed_blas.patch")
            patch.apply()
        os.chdir(cwd)







