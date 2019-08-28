.. highlight:: shell

============
Installation
============

Prerequisites
-------------

Create a virtual environment to avoid messing with python version.

Install conda: https://docs.conda.io/en/latest/miniconda.html

If you want to use sextractor to find sources you can create a python 3 environment. If you want to use pyraf you need to use python2.

For sextractor:

.. code-block:: console
 
    $ conda create -n gmadet python=3 numpy astropy astroquery matplotlib pandas 


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
    
    $ conda create -n gmadet27 python=2.7 iraf-all pyraf-all stsci



Activate the environment:

Sextractor:

.. code-block:: console
 
    $ conda activate gmadet 

Iraf:

.. code-block:: console
 
    $ conda activate iraf27

Install other libraries

.. code-block:: console
 
    $ pip install astroquery astroML requests h5py 


From sources
------------

The sources for gmadet can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/dcorre/gmadet

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/dcorre/gmadet/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/dcorre/gmadet
.. _tarball: https://github.com/dcorre/gmadet/tarball/master
