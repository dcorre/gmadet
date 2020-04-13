.. highlight:: shell

============
Installation
============

There are two ways for installing the project:
     * install all the dependencies yourself.
     * use a Docker to run the code.


The code is using some of the astromatic softwares that can be difficult to run on Windows. So, it is recommended to use the Docker for Windows users.


Prerequisites
-------------

For both ways you will need to get the project first. 

The sources for gmadet can be downloaded from the `Github repo`_.

* Either clone public repository:

  If git is not installed on your machine: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

.. code-block:: console

    $ git clone git://github.com/dcorre/gmadet

* Or download:

  Simply download this `tarball`_. Or through the console: 

.. code-block:: console

    $ curl  -OL https://github.com/dcorre/gmadet/tarball/master

Cloning the project allows to retrieve easily the future code updates. If you dowloaded the projet, you will need to download it again to retrieve future updates.


Installation with Docker
------------------------

The usage of a Docker allows to build an OS environment on your machine and thus avoid compatibility problems when running the code under Linux, Mac or Windows.   

* Install the Docker desktop for your OS: https://docs.docker.com/


Build the Docker image on your machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go to the gmadet/Docker/ directory, where the Dockerfile is.   
Type the following instruction:   

.. code-block:: console

   $ sudo docker image build -t gmadet:1.0 .
   
This will build a Ubuntu 18.04 image with the name gmadet with the version 1.0. This will install all the required libraries and softwares. This can take a few tens of minutes.

Run the Docker
^^^^^^^^^^^^^^

After building the image, you can run it in a terminal mode typing (from gmadet top directory):

.. code-block:: console

   $ sudo docker run -v $(pwd)/gmadet/:/home/newuser/gmadet/ --rm -it gmadet:1.0

This means that you run interactively in a bash terminal the Docker image. The -v option means that you mount a volume pointing to the gmadet directory on your host machine, so that you can exchange data between the Docker and your machine.

If the command $(pwd) does not work on your computer, just write the whole to the gmadet/gmadet/ repository on your computer.


Installation without Docker
---------------------------

Create a virtual environment to avoid messing with different python libraries version that could be already installed on your computer and required for other projects.

Install conda: https://docs.conda.io/en/latest/miniconda.html

If you want to use sextractor to find sources you can create a python 3 environment. If you want to use pyraf you need to use python2.


Python 3 environment with sextractor:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console
 
    $ conda create -n gmadet python=3 numpy scipy matplotlib astropy pandas shapely requests h5py sphinx sphinx_rtd_theme scikit-image


Python2 environment with Pyraf:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For pyraf, first install some 32bits libraries if your computer is a 64bits:

Debian >=7, Ubuntu >=14.04:

.. code-block:: console
 
    $ # If on Debian execute this first (not required on Ubuntu):
    $ sudo dpkg --add-architecture i386

    $ sudo apt-get update
    $ sudo apt-get install libc6:i386 libz1:i386 libncurses5:i386 libbz2-1.0:i386 libuuid1:i386 libxcb1:i386 libxmu6:i386

RHEL/CentOS >=6, Fedora >=14:

.. code-block:: console
 
    $ sudo yum install glibc.i686 zlib.i686 ncurses-libs.i686 bzip2-libs.i686 uuid.i686 libxcb.i686


.. code-block:: console
    
    $ conda create -n iraf27 python=2.7 iraf-all pyraf-all stsci shapely requests h5py sphinx sphinx_rtd_theme



Activate the environment:
^^^^^^^^^^^^^^^^^^^^^^^^^

Sextractor:

.. code-block:: console
 
    $ conda activate gmadet 

Iraf:

.. code-block:: console
 
    $ conda activate iraf27


Install other libraries
^^^^^^^^^^^^^^^^^^^^^^^

Once you have activated the environment, install the packages that are not available with conda using pip:

.. code-block:: console
 
    $ pip install lacosmic hjson voevent-parse xmltodict astroML regions photutils keras keras-vis
    $ pip install --pre astroquery

.. _Github repo: https://github.com/dcorre/gmadet
.. _tarball: https://github.com/dcorre/gmadet/tarball/master
