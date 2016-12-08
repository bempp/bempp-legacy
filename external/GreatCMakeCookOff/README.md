The Great CMake CookOff
=======================


This is a repository of usefull and less than usefull cmake recipes.  It is distributed under the
[MIT License](http://opensource.org/licenses/MIT)

Adding this repository to a cmake
=================================

The files in this repository can be added individually or as a whole to a project, as long as the
MIT copyright terms are followed. One possibility is to include this project as a [git
submodule](http://git-scm.com/docs/git-submodule).

However, the easiest method may well be to have this repository downloaded upon configuration of a
project. In that case, the file
[LookUp-GreatCMakeCookOff.cmake](https://github.com/UCL/GreatCMakeCookOff/tree/master/LookUp-GreatCMakeCookOff.cmake)
should be downloaded and inserted into the target project. It can then be included in the target
project's main `CMakeLists.txt` file:

```cmake
include(LookUp-GreatCMakeCookOff)
```

This will download the cook-off into the build directory right at configure time. Cook-off recipes
can then be used anywhere below that.

Another option is to point `CMake` towards the location on disk where a repo of the cook-off can be
found, or more explicitely, where the file `GreatCMakeCookOffConfig.cmake` can be found. This is
done with `cmake -DGreatCMakeCookOff_DIR=/path/to/cookoff/cmake ..`. Please note that this trick works
for any `CMake` project that defines `SomethingConfig.cmake` files.


Features
========

Please check the [wiki](https://github.com/UCL/GreatCMakeCookOff/wiki)
