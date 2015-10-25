Installation Instructions
=========================

Introduction
------------

`BEM++ <http://www.bempp.org>`_ is a C++/Python boundary element library.

Dependecies
-----------

To build and install BEM++ you will need to have the following installed:

*   Python 2 (2.? or newer) or 3. We recommend using `Enthought 
    Python <http://www.enthought.com/products/epd.php>`_ or `Anaconda
    Python <https://store.continuum.io/cshop/anaconda/>`_ as both come
    with MKL (for Anaconda, you will need to Accelerate package).

*   `Git version control system <http://git-scm.com/>`_.

*   `CMake <http://www.cmake.org/>`_.

Obtaining the Code
------------------

The code is available from https://github.com/bempp/bempp. You can
get the latest stable version be executing::

    git clone https://github.com/bempp/bempp.git

Compiling
---------

To compile BEM++, navigate to the directory where the source code has
downloaded::

    cd bempp

then make a folder called build and navigate into it::

    mkdir build
    cd build

BEM++ will compile in this directory. Next, make the configuration file
for the build with cmake::

    cmake ..

At this point, you can edit to locations of libraries by editing
CMakeCache.txt. Most importantly, make sure that the PYTHON_EXECUTABLE,
PYTHON_LIBRARY and PYTHON_INCLUDE_DIR variables are pointing to the
correct version of Python. To do this, you may find it easier to run
cmake within a python `virtual environment 
<https://virtualenv.pypa.io/en/latest/>`_.

Next, compile the BEM++ library::

    make

The installer will download any dependencies you do not have.

Running BEM++
-------------

To run a Python or IPython environment including BEM++, run::

    build/bin/bempp

or::

    build/bin/ibempp

If you run the following::

    export PATH="/path/to/bempp/build/bin:$PATH"

Then you will be able to run a Python or IPython environment including BEM++
by typing bempp or ibempp in terminal. You may find it useful to this to end of your .bashrc file.
