from os.path import basename, dirname, join, abspath
from setuptools import setup, Extension
from distutils.command.build import build as dBuild
from setuptools.command.install import install as dInstall
from setuptools.command.build_ext import build_ext as dBuildExt
from setuptools.command.bdist_egg import bdist_egg as dBuildDistEgg
from setuptools.command.sdist import sdist as dSDist
from setuptools.command.egg_info import egg_info as dEggInfo
from distutils.dir_util import mkpath

import sys

py_major = sys.version_info.major

source_dir = dirname(abspath(__file__))
package_dir = join(source_dir, 'pkg_install' + str(py_major))
mkpath(package_dir)

def cmake_cache_line(variable, value, type='STRING'):
    return "set(%s \"%s\" CACHE %s \"\")\n" % (variable, value, type)

def as_preload_file(name, info):
    """ Python information to cmake commandline """
    result = []
    if len(info.get('libraries', [])):
        libs = ('-l' + '-l'.join(info['libraries'])).rstrip().lstrip()
        if len(info.get('library_dirs', [])):
            libdirs =  ('-L' + '-L'.join(info['library_dirs'])).rstrip().lstrip()
        else: libdirs = ""
        result.append(cmake_cache_line("%s_LIBRARIES" % name,
            "%s %s" %(libdirs, libs)))
    if len(info.get('include_dirs', [])):
        incs = ';'.join(info['include_dirs']).rstrip().lstrip()
        result.append(cmake_cache_line("%s_INCLUDE_DIRS" % name, incs))
    return result

def cmake_executable():
    """ Path to cmake executable """
    from os.path import exists
    from os import environ
    from distutils.spawn import find_executable
    cmake = find_executable('cmake')
    if cmake is None and 'CASAPATH' in environ:
        # Tries to out-smart CASA.
        # Look places cmake might be that casa removes from path.
        directories = [
            join('/', 'usr', 'local', 'bin'),
            join(environ['HOME'], 'bin'),
            join(environ['HOME'], '.local', 'bin'),
            join(environ['HOME'], 'usr', 'bin'),
            join('/', 'sw', 'bin') # -- default Fink location
        ]
        for directory in directories:
            if exists(join(directory, 'cmake')):
                cmake = join(directory, 'cmake')
                break
    if cmake is None:
        raise RuntimeError('Could not find cmake executable in path')
    return cmake

class Build(dBuild):
    """ Build that runs cmake. """
    description = "Compiles BEM++ using cmake"
    user_options = dBuild.user_options + [
        ("external=", None, "Location for external packages")
    ]
    def initialize_options(self):
        self.external = None
        dBuild.initialize_options(self)

    def configure_cmdl(self, filename):
        """ Creates cmake command-line

            First puts variables into a cache file. This is safer that going through the
            command-line.
        """
        from sys import executable
        # other args
        other_args = [
            cmake_cache_line('PYTHON_EXECUTABLE', executable, 'PATH'),
            cmake_cache_line('NOEXPORT', 'TRUE', 'BOOL'),
            cmake_cache_line('PYPACKED', 'TRUE', 'BOOL'),
            cmake_cache_line('CMAKE_BUILD_TYPE', 'Release', 'STRING')
        ]
        if(self.external):
            other_args.extend([
                cmake_cache_line('EXTERNAL_ROOT', self.external, 'PATH'),
		cmake_cache_line('CMAKE_PREFIX_PATH', self.external, 'PATH')
            ])
        other_args.append('\n')

        with open(filename, 'w') as file: file.writelines(other_args)
        return ['-C%s' % filename]

    def _configure(self, build_dir):
        from distutils import log
        from distutils.spawn import spawn
        from os import chdir, getcwd

        current_dir = getcwd()
        mkpath(build_dir)
        command_line = self.configure_cmdl(join(build_dir, 'Variables.cmake'))
        log.info(
                "CMake: configuring with variables in %s "
                % join(build_dir, 'Variables.cmake')
        )
        cmake = cmake_executable()

        try:
            chdir(build_dir)
            spawn([cmake] + command_line + [source_dir])
        finally: chdir(current_dir)

    def _build(self, build_dir):
        from distutils import log
        from distutils.spawn import spawn
        from os import chdir, getcwd

        log.info("CMake: building in %s" % build_dir)
        current_dir = getcwd()
        cmake = cmake_executable()

        try:
            chdir(build_dir)
            spawn([cmake, '--build', '.'])
        finally: chdir(current_dir)

    def run(self):
        from os.path import abspath
	
        self.build_base = self.build_base + str(py_major)

        build_dir = join(dirname(abspath(__file__)), self.build_base)
        self._configure(build_dir)
        self._build(build_dir)
        self._install(build_dir, package_dir)
        try:
            prior = getattr(self.distribution, 'running_binary', False)
            dBuild.run(self)
        finally: self.distribution.running_binary = prior

    def _install(self, build_dir, install_dir):
        from distutils import log
        from distutils.sysconfig import PREFIX, get_python_lib
        from sys import version_info
        from os.path import abspath, relpath
        from os import chdir, getcwd

        libtopy = relpath(get_python_lib(), PREFIX)
        if len(libtopy) > 2 and libtopy[:1] == '..':
            libtopy = join(
                'lib',
                'python{0.major}.{0.minor}'.format(version_info),
                'site-packages'
            )

        current_cwd = getcwd()
        build_dir = abspath(build_dir)
        cmake = cmake_executable()
        install_dir = abspath(install_dir)
        log.info("CMake: Installing package to %s" % install_dir)
        try:
            chdir(build_dir)
            self.spawn([cmake,
                '-DPYTHON_PKG_DIR=\'%s\'' % install_dir,
                source_dir
            ])
            self.spawn([cmake, '--build', '.', '--target', 'install'])
        finally: chdir(current_cwd)
        self.distribution.running_binary = True

class Install(dInstall):
    def run(self):
        from distutils import log
        from os.path import abspath
        from os import chdir, getcwd
        self.distribution.run_command('build')
        current_cwd = getcwd()
        self.build_base = self.build_base + str(py_major)
        build_dir = join(dirname(abspath(__file__)), self.build_base)
        cmake = cmake_executable()
        pkg = abspath(self.install_lib)
        log.info("CMake: Installing package to %s" % pkg)
        try:
            chdir(build_dir)
            self.spawn([cmake,
                '-DPYTHON_PKG_DIR=\'%s\'' % pkg,
                '-DPYPACKED=TRUE',
                '..'
            ])
            self.spawn([cmake, '--build', '.', '--target', 'install'])
        finally: chdir(current_cwd)

        try:
            prior = getattr(self.distribution, 'running_binary', False)
            self.distribution.running_binary = True
            self.distribution.have_run['egg_info'] = 0
            dInstall.run(self)
        finally: self.distribution.running_binary = prior

class BuildExt(dBuildExt):
    def __init__(self, *args, **kwargs):
        dBuildExt.__init__(self, *args, **kwargs)
    def run(self):pass

class BuildDistEgg(dBuildDistEgg):
    def __init__(self, *args, **kwargs):
        dBuildDistEgg.__init__(self, *args, **kwargs)
    def run(self):

        try:
            prior = getattr(self.distribution, 'running_binary', False)
            self.distribution.running_binary = True
            self.run_command('build')
            dBuildDistEgg.run(self)
        finally: self.distribution.running_binary = prior

class EggInfo(dEggInfo):
    def __init__(self, *args, **kwargs):
        dEggInfo.__init__(self, *args, **kwargs)
    def run(self):
        from setuptools.command.egg_info import manifest_maker
        from os import listdir
        which_template = 'MANIFEST.source.in'

        dist = self.distribution
        old_values = dist.ext_modules, dist.ext_package, \
            dist.packages, dist.package_dir
        if len(listdir(package_dir)) != 0  \
            and getattr(self.distribution, 'running_binary', False):
            which_template = 'MANIFEST.binary.in'
        else:
            dist.ext_modules, dist.ext_package = None, None
            dist.packages, dist.package_dir = None, None

        try:
            old_template = manifest_maker.template
            manifest_maker.template = which_template
            dEggInfo.run(self)
        finally:
            manifest_maker.template = old_template
            dist.ext_modules, dist.ext_package = old_values[:2]
            dist.packages, dist.package_dir = old_values[2:]

class SDist(dSDist):
    def __init__(self, *args, **kwargs):
        dSDist.__init__(self, *args, **kwargs)
    def run(self):
        dist = self.distribution
        try:
            old_values = dist.ext_modules, dist.ext_package, \
                dist.packages, dist.package_dir
            dist.ext_modules, dist.ext_package = None, None
            dist.packages, dist.package_dir = None, None
            dSDist.run(self)
        finally:
            dist.ext_modules, dist.ext_package = old_values[:2]
            dist.packages, dist.package_dir = old_values[2:]
setup(
    name = "bempp",
    version = "3.1.5",

    setup_requires = ['numpy', 'scipy', 'cython>=0.23', 'mpi4py'],
    install_requires = ['numpy', 'scipy', 'cython>=0.23', 'mpi4py'],
    platforms = ['GNU/Linux','Unix','Mac OS-X'],

    zip_safe = False,
    cmdclass = {
        'build': Build, 'install': Install,
        'build_ext': BuildExt, 'bdist_egg': BuildDistEgg,
        'egg_info': EggInfo
    },

    author = "Timo Betcke",
    author_email = "timo.betcke@gmail.com",
    description = "The BEM++ boundary element library",
    license = "MIT",
    url = "https://bitbucket.com/bemppsolutions/bempp",
    ext_modules = [Extension('bempp.core', [])],
    ext_package = 'bempp',
    packages = ['bempp'],
    package_dir = {
        'bempp': join(basename(package_dir), 'bempp')
    },
    include_package_data=True,

    keywords= "boundary element method",
    classifiers = [
         'Development Status :: 5 - Production/Stable',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: MIT License',
         'Operating System :: MacOS X',
         'Operating System :: Linux',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3.4',
         'Topic :: Scientific/Engineering',
         'Topic :: Scientific/Engineering :: Mathematics',
    ])
    #long_description = open(join(dirname(__file__), 'README'), 'r').read()

