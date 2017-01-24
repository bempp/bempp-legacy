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

def split_regex(lines, pattern):
    """ Split lines into matches and lines that do not match """
    haves, have_nots = [], []
    for line in lines:
        m = pattern.match(line)
        if m is not None: haves.append(m)
        else: have_nots.append(line)
    return haves, have_nots

def parse_ldd_output(output):
    """Search ldd output for MKL dependencies.

    Return (mkl_dirs, mkl_libs)."""
    from os.path import dirname, basename
    from re import compile

    # like libvtkWidgets.so.pv3.12 => not found
    re1 = compile(r"\s*(.+) => not found")
    # like "libz.so.1 => /usr/lib/libz.so.1 (0xb7d60000)"
    re2 = compile(r"\s*(.+) => (.+) \(.+\)")
    # like /lib/ld-linux.so.2 (0x80000000)
    re3 = compile(r"\s*(/.+) \(.+\)")
    # # like linux-gate.so.1 =>  (0xffffe000)
    # (don't know what to do with this form)
    # re4 = re.compile(r"\s*(/.+) => \(.+\)")

    re_fname = compile(r"lib(mkl.*|iomp.*)\.(so|dylib)(\.[^ ]*)?")
    mkl_dirs = []
    mkl_libs = []
    re1_match, output_lines = split_regex(output.splitlines(), re1)
    for m in re1_match:
        fname = m.group(1)
        m_fname = re_fname.match(fname)
        if m_fname:
            # can't do better than this since the full path is unknown...
            mkl_libs.append("-l"+m_fname.group(1))
            print ("Warning: NumPy MKL dependency '"+fname+"' not found")

    re2_match, output_lines = split_regex(output_lines, re2)
    for m in re2_match:
        fname = m.group(1)
        path = m.group(2)
        m_fname = re_fname.match(fname)
        if m_fname:
            mkl_libs.append(path)
            mkl_dirs.append(dirname(path))

    for m in split_regex(output_lines, re3)[0]:
        path = m.group(1)
        fname = basename(path)
        m_fname = re_fname.match(fname)
        if m_fname:
            mkl_libs.append(path)
            mkl_dirs.append(dirname(path))

    return set(mkl_dirs), set(mkl_libs)

def parse_otool_output(output):
    """Search otool output for MKL dependencies.

    Return (mkl_dirs, mkl_libs)."""
    from re import compile
    from os.path import join, dirname, split as split_path, abspath, basename
    import numpy
    import sys

    # like "@rpath/libmkl_intel.dylib (compatibility version 0.0.0, current version 0.0.0)"
    re1 = compile(r"\s*@rpath/lib/(.+) \(.+\)")
    # like "@loader_path/libmkl_intel.dylib (compatibility version 0.0.0, current version 0.0.0)"
    re2 = compile(r"\s*@loader_path/(.+) \(.+\)")
    # like "/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 111.0.0)"
    re3 = compile(r"\s*(.+) \(.+\)")
    re_fname = compile(r"lib(mkl.*|iomp.*)\.(so|dylib)(\.[^ ]*)?")
    # we assume for now that @rpath == <sys.prefix>/lib
    prefix_dir = getattr(sys, 'base_prefix', sys.prefix)
    sys_lib_dir = join(prefix_dir, "lib")

    mkl_dirs, mkl_libs = [], []
    re1_match, output_lines = split_regex(output.splitlines(), re1)
    for m in re1_match:
        fname = m.group(1)
        m_fname = re_fname.match(fname)
        if m_fname:
            # we assume that @rpath is equal to sys.prefix
            mkl_libs.append(join(sys_lib_dir, m.group(1)))
            mkl_dirs.append(sys_lib_dir)

    re2_match, output_lines = split_regex(output_lines, re2)
    for m in re2_match:
        full_path = join(dirname(numpy.__file__), 'linalg', m.group(1))
        fpath, fname = split_path(abspath(full_path))
        m_fname = re_fname.match(fname)
        if m_fname:
            mkl_libs.append(full_path)
            mkl_dirs.append(fpath)

    for m in split_regex(output_lines, re3)[0]:
        path = m.group(1)
        fname = basename(path)
        m_fname = re_fname.match(fname)
        if m_fname:
            mkl_libs.append(path)
            mkl_dirs.append(dirname(path))

    return set(mkl_dirs), set(mkl_libs)

def get_mkl_dirs_and_libs_like_numpy():
    from os.path import dirname, join, isfile
    from sys import platform
    from os import environ
    from subprocess import check_output
    from glob import glob
    try:
        import numpy
    except ImportError:
        raise Exception("MKL autodetection failed: NumPy could not "
                        "be imported. Specify MKL location manually")

    lapack_lite_path = glob(join(dirname(numpy.__file__), 'linalg', 'lapack_lite*.so'))[0]
    if not isfile(lapack_lite_path):
        raise Exception("MKL autodetection failed: '"+lapack_lite_path+
                        "' is not a file. Specify MKL location manually")
    if platform.startswith('darwin'):
        otool_output = check_output(['otool','-L',lapack_lite_path]).decode("utf-8")
        mkl_dirs,mkl_libs = parse_otool_output(otool_output)
    else: # 'linux' -- we've checked that its 'darwin' or 'linux' before
        env = environ.copy()
        # Fix for Canopy 1.0.3: ensure that LD_LIBRARY_PATH contains the
        # appdata/canopy-*/lib directory
        if "VENV_LD_LIBRARY_PATH" in env:
            env["LD_LIBRARY_PATH"] = env["VENV_LD_LIBRARY_PATH"]
        ldd_output = check_output(['ldd',lapack_lite_path],env=env).decode("utf-8")
        mkl_dirs,mkl_libs = parse_ldd_output(ldd_output)
    return set(mkl_dirs), set(mkl_libs)

if __name__ == '__main__':
    # Output that CMake can parse, or error
    from sys import exit
    directories, libraries = get_mkl_dirs_and_libs_like_numpy()
    if not len(libraries): exit(1)
    result = ";".join(libraries)
    print(result)
