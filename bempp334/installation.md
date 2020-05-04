---
title: Installing Bempp 3.3.4
---

On this page, you can find information about installing Bempp 3.3.4.
This is the old, C++ based version of the library.
You can find information about installing the newest version of Bempp ({{ site.bemppversion }}) [here](../installation.md).


## Running with Docker
The easiest way to get started using Bempp is to use the Bempp Docker image.
This provides a Jupyter notebook server with a full Bempp environment, capable of running Python 2 or Python 3 notebooks.
For Python 2 notebooks, FEniCS is provided to solve FEM/BEM coupled problems.
The container is based on Ubuntu 17.04 Zesty.

Once you have installed Docker the Bempp image can be pulled using

```bash
docker pull bempp/notebook
```

This command will download the latest image to your machine.
If a new version of Bempp is released just repeat this command to obtain the latest version.
Depending on the installation (such as on a typical Linux system) you may need to run the above command as root user. To do this simply run

```bash
sudo docker pull bempp/notebook
```

The Bempp notebook server can then be started using the following command:

```bash
docker run -it --rm -v $(pwd):/home/bempp/work -p 8888:8888 bempp/notebook
```

Again, prepend this command with `sudo` if necessary.
This command will automatically start the notebook server on localhost:8888.
Simply enter this address in your web browser and a Jupyter Notebook server with Bempp appears.

### Password protection
By default, the notebook server accepts any connection on port 8888.
If you are in a multi-user environment, or desire additional security for any other reason, you can use the `start-notebook.sh` script to provide options to the notebook server.

For example, to start a password protected Notebook server, first generate a password inside an IPython session with the following commands:

```python
from notebook.auth import passwd
passwd()
```

This will ask you for a password, then give you its SHA-1 hash in the form `'sha1:...'`. We can then start the notebook server with the following command (prepend `sudo` if necessary):

```bash
docker run -it --rm -v $(pwd):/home/bempp/work -p 8888:8888 bempp/notebook start-notebook.sh --NotebookApp.password='sha1:...'
```

The page at localhost:8888 will now prompt you for the password before allowing you to run notebooks.

### Visualization in notebooks

Please note that in order to use interactive visualization in this setup you need to enable the IPython viewer after importing bempp into a Jupyter notebook, such as by

```python
import bempp.api
bempp.api.set_ipython_notebook_viewer()
```

This uses the WebGL feature of modern web browsers to enable the plot methods for grids and grid functions within the Jupyter notebook.
This solution is meant only for prototyping.
For large meshes and advanced visualization options it is recommended to export grids or grid functions as Gmsh files and to view them externally.

## Building from source
You can find the source code of Bempp 3.4.4 on [Bitbucket](https://bitbucket.org/bemppsolutions/bempp/src/master/).

Once you have downloaded the source code, you can build it by running the following commands in terminal:

```bash
mkdir build
cd build
cmake ..
make
```
