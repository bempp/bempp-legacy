##########################

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

##########################
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

##########################
def writeOptions(root,config):
    """Write out options in a format suitable for a shell script"""
   
    fname=root+"/options.cfg"
    f=open(fname,'w')
    f.write("#This file is created automatically by bemppInstall.py\n")
    for section in config.sections():
        for option in config.items(section):
            f.write(section+"_"+option[0]+"="+"\""+option[1]+"\""+"\n")
    f.close()

############################

def setDefaultConfigOption(config,section,option,value):
    """Enter a default option into the ConfigParser object 'config'. If option already exists returns
       the existing value, otherwise the value 'value'.
    """

    if config.has_option(section,option):
        return config.get(section,option)
    else:
        if not config.has_section(section): config.add_section(section)
        config.set(section,option,value)
        return value

