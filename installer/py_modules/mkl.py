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

import glob,os,re,shutil,subprocess,sys,urllib
from py_modules import tools

def find_file_in_dirs(fname, dirs):
    for d in dirs:
        path = os.path.join(d, fname)
        if os.path.isfile(path):
            return path
    return None

def parse_ldd_output(output):
    """Search ldd output for MKL dependencies.

    Return (mkl_dirs, mkl_libs)."""

    # like libvtkWidgets.so.pv3.12 => not found
    re1 = re.compile(r"\s*(.+) => not found")
    # like "libz.so.1 => /usr/lib/libz.so.1 (0xb7d60000)"
    re2 = re.compile(r"\s*(.+) => (.+) \(.+\)")
    # like /lib/ld-linux.so.2 (0x80000000)
    re3 = re.compile(r"\s*(/.+) \(.+\)")
    # # like linux-gate.so.1 =>  (0xffffe000)
    # (don't know what to do with this form)
    # re4 = re.compile(r"\s*(/.+) => \(.+\)")

    re_fname = re.compile(r"lib(mkl.*|iomp.*)\.(so|dylib)(\.[^ ]*)?")
    mkl_dirs = []
    mkl_libs = []
    for l in output.splitlines():
        m = re1.match(l)
        if m:
            fname = m.group(1)
            m_fname = re_fname.match(fname)
            if m_fname:
                # can't do better than this since the full path is unknown...
                mkl_libs.append("-l"+re_fname.group(1))
                print "Warning: NumPy MKL dependency '"+fname+"' not found"
            continue
        m = re2.match(l)
        if m:
            fname = m.group(1)
            path = m.group(2)
            m_fname = re_fname.match(fname)
            if m_fname:
                mkl_libs.append(path)
                mkl_dirs.append(os.path.dirname(path))
            continue
        m = re3.match(l)
        if m:
            path = m.group(1)
            fname = os.path.basename(path)
            m_fname = re_fname.match(fname)
            if m_fname:
                mkl_libs.append(path)
                mkl_dirs.append(os.path.dirname(path))
    return mkl_dirs,mkl_libs

def parse_otool_output(output):
    """Search otool output for MKL dependencies.

    Return (mkl_dirs, mkl_libs)."""

    # like "@rpath/libmkl_intel.dylib (compatibility version 0.0.0, current version 0.0.0)"
    re1 = re.compile(r"\s*@rpath/(.+) \(.+\)")
    # like "/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 111.0.0)"
    re2 = re.compile(r"\s*(.+) \(.+\)")

    re_fname = re.compile(r"lib(mkl.*|iomp.*)\.(so|dylib)(\.[^ ]*)?")
    # we assume for now that @rpath == <sys.prefix>/lib
    if hasattr(sys,'base_prefix'):
        prefix_dir = sys.base_prefix
    else:
        prefix_dir = sys.prefix
    sys_lib_dir = os.path.join(prefix_dir, "lib") 
    mkl_dirs = []
    mkl_libs = []
    for l in output.splitlines():
        m = re1.match(l)
        if m:
            fname = m.group(1)
            m_fname = re_fname.match(fname)
            if m_fname:
                # we assume that @rpath is equal to sys.prefix
                mkl_libs.append(os.path.join(sys_lib_dir,m.group(1)))
                mkl_dirs.append(sys_lib_dir)
            continue
        m = re2.match(l)
        if m:
            path = m.group(1)
            fname = os.path.basename(path)
            m_fname = re_fname.match(fname)
            if m_fname:
                mkl_libs.append(path)
                mkl_dirs.append(os.path.dirname(path))
    return mkl_dirs,mkl_libs

def get_mkl_dirs_and_libs_like_numpy(config, lib_dir, extension):
    try:
        import numpy
    except ImportError:
        raise Exception("MKL autodetection failed: NumPy could not "
                        "be imported. Specify MKL location manually")
    numpy_path = numpy.__file__ # path to numpy/__init__.pyc
    numpy_dir = os.path.dirname(numpy_path)
    lapack_lite_path = os.path.join(numpy_dir,
                                    "linalg/lapack_lite.so")
    if not os.path.isfile(lapack_lite_path):
        raise Exception("MKL autodetection failed: '"+lapack_lite_path+
                        "' is not a file. Specify MKL location manually")
    if sys.platform.startswith('darwin'):
        otool_output = tools.check_output(['otool','-L',lapack_lite_path])
        mkl_dirs,mkl_libs = parse_otool_output(otool_output)
    else: # 'linux' -- we've checked that its 'darwin' or 'linux' before
        ldd_output = tools.check_output(['ldd',lapack_lite_path])
        mkl_dirs,mkl_libs = parse_ldd_output(ldd_output)
    return mkl_dirs,mkl_libs

def get_mkl_dirs_and_libs_installed(config, lib_dir, extension):
    mkl_dirs = config.get('MKL','dir')
    mkl_dirs = mkl_dirs.split(";")
    mkl_dirs = [s.strip() for s in mkl_dirs]
    if not mkl_dirs:
        raise Exception("Option 'dir' in section 'MKL' "
                        "must not be empty")
    for d in mkl_dirs:
        if not os.path.isdir(d):
            raise Exception("Option 'dir' in section 'MKL' "
                            "is invalid: '"+d+"' is not a directory")
    mkl_libs = config.get('MKL','lib').split(";")
    mkl_libs = [s.strip() for s in mkl_libs]
    return mkl_dirs,mkl_libs

def get_mkl_dirs_and_libs_redistributable(config, lib_dir, extension):
    mkl_tarball = config.get('MKL','tarball')
    mkl_tarball = tools.normalizePath(config,mkl_tarball)
    print 'Extracting MKL redistributables'
    tools.extract_file(mkl_tarball,lib_dir)
    mkl_dirs = [lib_dir]
    mkl_libs = config.get('MKL','lib').split(";")
    mkl_libs = [s.strip() for s in mkl_libs]
    return mkl_dirs,mkl_libs

# Find MKL library files and create symbolic links to them
# in the lib_dir directory
def create_symlinks(lib_dir, extension, mkl_dirs, mkl_libs):
    # First, remove any old symlinks to known MKL libraries
    # from the prefix directory
    for path in glob.glob(os.path.join(lib_dir,"libmkl")+"*"+extension):
        os.unlink(path)
    path = os.path.join(lib_dir,"libiomp5"+extension)
    if os.path.exists(path):
        os.unlink(path)

    # Create symlinks to all libs listed in mkl_libs
    for l in mkl_libs:
        if l.startswith("-l"):
            fname = "lib"+l[2:]+extension
            path = find_file_in_dirs(fname, mkl_dirs)
            if not path:
                raise Exception("MKL library '"+fname+"' not found in any"
                                "of these directories: '"+
                                "', '".join(mkl_dirs)+"'")
        elif os.path.isfile(l):
            fname = os.path.basename(l)
            path = l
        else:
            raise Exception("'"+l+"' is neither a file nor "
                            "a linker '-l' option")
        new_path = os.path.join(lib_dir,fname)
        if os.path.exists(new_path):
            os.unlink(new_path) # maybe it's a leftover symlink from a previous
                                # run, with a nonstandard name
        os.symlink(path,new_path)
    # Symlink to all other libmkl* files found in mkl_dirs
    for d in mkl_dirs:
        for path in glob.glob(os.path.join(d,"libmkl")+"*"+extension):
            fname = os.path.basename(path)
            new_path = os.path.join(lib_dir,fname)
            if not os.path.exists(new_path):
                os.symlink(path,new_path)
    # Symlink to the first libiomp5* file found in mkl_dirs
    fname = "libiomp5"+extension
    for d in mkl_dirs:
        path = os.path.join(d,fname)
        if os.path.exists(path):
            new_path = os.path.join(lib_dir,fname)
            if not os.path.exists(new_path):
                os.symlink(path,new_path)
            break

def get_linker_args(lib_dir, extension, mkl_dirs, mkl_libs):
    linker_args = []
    for l in mkl_libs:
        if l.startswith("-l"):
            path = os.path.join(lib_dir,"lib"+l[2:]+extension)
            if os.path.isfile(path):
                linker_args.append(path)
            else:
                linker_args.append(l)
        else:
            path = os.path.join(lib_dir,os.path.basename(l))
            linker_args.append(path)
    return linker_args

def download(root,config):
    pass

def prepare(root,config):
    prefix = config.get('Main','prefix')
    lib_dir = prefix+"/bempp/lib"

    if sys.platform.startswith('darwin'):
        extension = ".dylib"
    elif sys.platform.startswith('linux'):
        extension = ".so"
    else:
        raise Exception("Unsupported platform: '"+sys.platform+"'")

    enable_mkl = tools.to_bool(tools.setDefaultConfigOption(config,'MKL',
                                                            'enable_mkl','no'))
    if enable_mkl:
        if config.has_option('BLAS','lib') and config.get('BLAS','lib'):
            print ("Warning: contents of the option 'lib' in section 'BLAS' "
                   "will be ignored, since enable_mkl is set")
        if config.has_option('LAPACK','lib') and config.get('LAPACK','lib'):
            print ("Warning: contents of the option 'lib' in section 'LAPACK' "
                   "will be ignored, since enable_mkl is set")

        # Variables:
        # - mkl_dirs: list of directories with all MKL libraries
        # - mkl_libs: list of MKL libraries to link against, either as
        #   full paths or as '-l...' linker commands
        # - mkl_linker_args: list of MKL libraries to link against,
        #   either as full paths to files contained in lib_dir or as
        #   '-l...' linker commands

        # If mkl_source is not set, we'll get a sensible exception message
        mkl_source = config.get('MKL','source').lower()
        if mkl_source == "like_numpy":
            mkl_dirs,mkl_libs = get_mkl_dirs_and_libs_like_numpy(
                config,lib_dir,extension)
            if not mkl_libs:
                raise Exception(
                    "Your NumPy package does not seem to be linked to MKL. "
                    "Are you sure\nyou are using the correct version of Python? "
                    "If so, set the 'source'\noption in the 'MKL' section of "
                    "your configuration file to 'installed'\nand use the "
                    "options 'dir' and 'lib' to specify manually the location "
                    "of your\nMKL libraries.")
            create_symlinks(lib_dir,extension,mkl_dirs,mkl_libs)
        elif mkl_source == 'installed':
            mkl_dirs,mkl_libs = get_mkl_dirs_and_libs_installed(
                config,lib_dir,extension)
            create_symlinks(lib_dir,extension,mkl_dirs,mkl_libs)
        elif mkl_source == 'redistributable':
            mkl_dirs,mkl_libs = get_mkl_dirs_and_libs_redistributable(
                config,lib_dir,extension)
            # no need to create symlinks
        else:
            raise Exception("Option 'mkl_source' in section 'MKL' must be "
                            "either 'installed', 'redistributable' or "
                            "'like_numpy'")

        mkl_linker_args = get_linker_args(lib_dir,extension,mkl_dirs,mkl_libs)
        blas_libs = ";-lpthread;".join(mkl_linker_args)+";-lpthread"
        lapack_libs = blas_libs

        tools.setDefaultConfigOption(config,'BLAS','lib',blas_libs,
                                     overwrite=True)
        tools.setDefaultConfigOption(config,'LAPACK','lib',lapack_libs,
                                     overwrite=True)
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









