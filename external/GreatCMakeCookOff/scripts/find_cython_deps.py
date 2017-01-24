from re import compile

from_extern_re = "^\s*cdef\s+extern\s+from\s+"
from_extern_re += "(?P<quote>'|\")(?P<filename>\S+)(?P=quote)\s*:"
from_extern_re = compile(from_extern_re)

from_cimport_re = compile("^\s*from\s+(?P<filename>\S+)\s+cimport")


def check_c_includes(filename, includes):
    """ Check whether file exist in include dirs """
    from os.path import exists, isfile, join
    for directory in includes:
        path = join(directory, filename)
        if exists(path) and isfile(path):
            return path


def check_cython_includes(filename, includes):
    """ Check whether file exist in include dirs """
    from os.path import exists, isfile, join
    for directory in includes:
        path = join(directory, filename) + ".pxd"
        if exists(path) and isfile(path):
            return path
        path = join(directory, *filename.split('.')) + ".pxd"
        if exists(path) and isfile(path):
            return path


def check_extern_dep(line, results, includes):
    from os.path import abspath
    found = from_extern_re.match(line)
    if found is not None:
        actual = check_c_includes(found.group('filename'), includes)
        if actual is not None:
            results.add(abspath(actual))


def check_cimport_dep(line, results, includes):
    from os.path import abspath
    found = from_cimport_re.match(line)
    if found is not None:
        actual = check_cython_includes(found.group('filename'), includes)
        if actual is not None and actual not in results:
            results.add(abspath(actual))
            cython_deps(actual, includes, results)


def cython_deps(path, includes, results=None):
    from os.path import splitext, split, exists, isfile, join, abspath
    if results is None:
        results = set([abspath(path)])

    path = abspath(path)
    directory, filename = split(path)
    filename, extension = splitext(filename)
    if extension == '.pyx':
        newpath = join(directory, filename + '.pxd')
        if exists(newpath) and isfile(newpath):
            results.add(abspath(newpath))
            cython_deps(newpath, includes, results)

    # Generated files may not exist
    if not exists(path): 
        return results

    with open(path, 'r') as file:
        for line in file:
            check_extern_dep(line, results, includes)
            check_cimport_dep(line, results, includes)
    return results

if __name__ == '__main__':
    from sys import argv
    from os.path import abspath
    includes = set([abspath(u) for u in argv[2:]])
    deps = cython_deps(argv[1], includes)
    print(";".join(deps))
