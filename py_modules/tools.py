import os,subprocess,sys, shutil
from subprocess import CalledProcessError


# The following function extracts the tar gz files. It is taken from
# http://code.activestate.com/recipes/576714-extract-a-compressed-file/

def extract_file(path, to_directory='.'):
    import tarfile,zipfile,os
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

# Taken from 
# http://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python

def to_bool(value):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("on", "yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("off", "no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))

def writeOptions(root,config):
    """Write out options in a format suitable for a shell script"""
   
    fname=root+"/.options.cfg"
    f=open(fname,'w')
    f.write("#This file is created automatically by bemppInstall.py\n")
    for section in config.sections():
        for option in config.items(section):
            f.write(section+"_"+option[0]+"="+"\""+option[1]+"\""+"\n")
    f.close()


def setDefaultConfigOption(config,section,option,value, overwrite=False):
    """Enter a default option into the ConfigParser object 'config'. If option already exists returns
       the existing value, otherwise the value 'value'.
    """

    if config.has_option(section,option) and not overwrite:
        return config.get(section,option)
    else:
        if not config.has_section(section): config.add_section(section)
        config.set(section,option,value)
        return value

def pythonInfo():
    """Return a tuple (exe,lib,include) with the paths of the Python Interpeter, Python library and include directory"""

    import sys,os
    
    exe = sys.executable
    lib_no_suffix = sys.prefix+"/lib/libpython"+str(sys.version_info.major)+"."+str(sys.version_info.minor)
    if sys.platform.startswith('darwin'):
        lib = lib_no_suffix+".dylib"
    elif sys.platform.startswith('linux'):
        lib = lib_no_suffix+".so"
    else:
        raise Exception("Platform not supported")
    if not os.path.isfile(lib):
        lib = lib_no_suffix+".a"
        if not os.path.isfile(lib):
            raise Exception("Could not find Python library in "+sys.prefix+"/lib/")
    include = sys.prefix+"/include/python"+str(sys.version_info.major)+"."+str(sys.version_info.minor)
    return (exe,lib,include)

def download(fname,url,dir):
    """Download a file from a url into the directory dir if the file does not already exist"""
    from py_modules.urlgrabber import urlgrab,progress
    if not os.path.isfile(dir+"/"+fname):
        urlgrab(url,filename=dir+"/"+fname,progress_obj=progress.TextMeter(),reget='simple')

def checkCreateDir(dir):
    """Create a directory if it does not yet exist"""
    if not os.path.isdir(dir):
        os.mkdir(dir)

def to_int(s,min_val=1):
    """Convert string s to integer. if int(s)<min_val return min_val"""
    i = int(s)
    if i>min_val:
        return i
    else:
        return min_val

def testBlas(root,config):
    """Test if BLAS functions correctly and set whether G77 calling convention is needed"""
    cwd = os.getcwd()
    blas = config.get('BLAS','lib')
    cxx=config.get('Main','cxx')
    cc=config.get('Main','cc')
    cflags = config.get('Main','cflags')
    cxxflags = config.get('Main','cxxflags')
    prefix = config.get('Main','prefix')
    fnull = open(os.devnull,'w')

    if sys.platform.startswith('darwin'):
        ld_path = "export DYLD_LIBRARY_PATH="+prefix+"/bempp/lib:$DYLD_LIBRARY_PATH; "
    elif sys.platform.startswith('linux'):
        ld_path = "export LD_LIBRARY_PATH="+prefix+"/bempp/lib:$LD_LIBRARY_PATH; "
    else:
        raise Exception("Wrong architecture.")


    os.mkdir(root+"/test_blas/build")
    os.chdir(root+"/test_blas/build")
    config_string = "CC="+cc+" CXX="+cxx+" CFLAGS='"+cflags+"' CXXFLAGS='"+cxxflags+"' cmake -D BLAS_LIBRARIES:STRING=\""+blas+"\" .."
    try:
        subprocess.check_call(config_string,shell=True,stdout=fnull,stderr=fnull)
        subprocess.check_call("make",shell=True,stdout=fnull,stderr=fnull)
    except:
        fnull.close()
        raise Exception("Building BLAS tests failed. Please check your compiler and BLAS library settings")
    try:
        subprocess.check_call(ld_path+ "./test_blas",shell=True,stdout=fnull,stderr=fnull)
    except:
        os.chdir(cwd)
        fnull.close()
        raise Exception("BLAS is not working correctly. Please check your libraries and your DYLD_LIBRARY_PATH or LD_LIBRARY_PATH settings.")
    # Test for zdotc convention
    g77 = False
    try:
        subprocess.check_call(ld_path+"./test_zdotc",shell=True,stdout=fnull,stderr=fnull)
    except:
        g77 = True
    if g77:
        try:
            subprocess.check_call(ld_path+"./test_zdotc_g77",shell=True,stdout=fnull,stderr=fnull)
        except:
            os.chdir(cwd)
            fnull.close()
            raise Exception("ZDOTC works neither in G77 nor in GFortran modus. Please check your BLAS libraries")
        setDefaultConfigOption(config,'AHMED','with_g77','ON',overwrite=True)
    else:
        setDefaultConfigOption(config,'AHMED','with_g77','OFF',overwrite=True)
    os.chdir(cwd)
    print "BLAS configuration successfully completed."
    
def testLapack(root,config):
    """Test if BLAS functions correctly and set whether G77 calling convention is needed"""
    cwd = os.getcwd()
    blas = config.get('BLAS','lib')
    lapack = config.get('LAPACK','lib')
    cxx=config.get('Main','cxx')
    cc=config.get('Main','cc')
    cflags = config.get('Main','cflags')
    cxxflags = config.get('Main','cxxflags')
    prefix = config.get('Main','prefix')
    fnull = open(os.devnull,'w')

    if sys.platform.startswith('darwin'):
        ld_path = "export DYLD_LIBRARY_PATH="+prefix+"/bempp/lib:$DYLD_LIBRARY_PATH; "
    elif sys.platform.startswith('linux'):
        ld_path = "export LD_LIBRARY_PATH="+prefix+"/bempp/lib:$LD_LIBRARY_PATH; "
    else:
        raise Exception("Wrong architecture.")


    os.mkdir(root+"/test_lapack/build")
    os.chdir(root+"/test_lapack/build")
    config_string = "CC="+cc+" CXX="+cxx+" CFLAGS='"+cflags+"' CXXFLAGS='"+cxxflags+"' cmake -D BLAS_LIBRARIES:STRING=\""+blas+"\" -D LAPACK_LIBRARIES=\""+lapack+"\" .."
    try:
        subprocess.check_call(config_string,shell=True,stdout=fnull,stderr=fnull)
        subprocess.check_call("make",shell=True,stdout=fnull,stderr=fnull)
    except:
        fnull.close()
        raise Exception("Building LAPACK tests failed. Please check your compiler and BLAS library settings")
    try:
        subprocess.check_call(ld_path+ "./test_lapack",shell=True,stdout=fnull,stderr=fnull)
    except:
        os.chdir(cwd)
        fnull.close()
        raise Exception("LAPACK is not working correctly. Please check your libraries and your DYLD_LIBRARY_PATH or LD_LIBRARY_PATH settings.")
    os.chdir(cwd)
    print "LAPACK configuration successfully completed."

def checkDeleteDirectory(dir):
    """Delete directory dir if it exists"""
    if os.path.isdir(dir): shutil.rmtree(dir)

def checkDeleteFile(s):
    """Delete file given by path in string s if it does not exist"""
    if os.path.exists(s):
        os.remove(s)

def cleanUp(root,config):
    """Clean up so that a rerun of the installer is possible"""

    prefix=config.get('Main','prefix')
    checkDeleteDirectory(prefix+"/bempp")
    checkDeleteDirectory(root+"/test_blas/build")
    checkDeleteDirectory(root+"/test_lapack/build")

    
        
    
    

