.. highlight:: shell

============
Installation
============

There are two ways for installing the project:
     * install all the dependencies yourself.
     * use a Docker image to run the code.

The Docker image contains all the C dependencies and python packages, that can be sometimes very time consuming to install depending on your OS. So I recommend using the Docker image.
The code is using some of the astromatic softwares that can be difficult to run on Windows. So, it is strongly advised  to use the Docker image for Windows users. 


Prerequisites
-------------

For both ways you will need to get the project first. 

The sources for gmadet can be downloaded from the `Github repo`_.

* Either clone public repository:

  If git is not installed on your machine, see: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

  Then clone the project:

.. code-block:: console

    $ git clone git://github.com/dcorre/gmadet

* Or download:

  Simply download this `tarball`_. Or through the console: 

.. code-block:: console

    $ curl  -OL https://github.com/dcorre/gmadet/tarball/master

Cloning the project allows to retrieve easily future code updates. If you downloaded the projet, you will need to download it again to retrieve future updates.


Installation with Docker
------------------------

The usage of a Docker allows to build an OS environment on your machine and thus avoid compatibility problems when running the code under Linux, Mac or Windows. If you have not Docker installed on your machine install it first.   

* Install the Docker desktop for your OS: https://docs.docker.com/get-docker/

* To run Docker without appending sudo, type:

.. code-block:: console
   
   $ sudo groupadd docker
   $ sudo usermod -aG docker $USER

Log out and log back in so that your group membership is re-evaluated. For more information see https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user.

You can test that Docker is installed correctly and can be run without sudo:

.. code-block:: console

   $ docker run hello-world


Download the gmadet Docker image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To retrieve the Docker image:

.. code-block:: console

   $ docker pull dcorre/gmadet

Check that it appears in the list of images:

.. code-block:: console

   $ docker images


Installation without Docker
---------------------------

I advise to create a virtual environment to avoid messing with different python libraries version that could be already installed on your computer and required for other projects.

Install conda: https://docs.conda.io/en/latest/miniconda.html

You can also install everything with pip if you prefer not to use conda.

Python 3 environment:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console
 
    $ conda create -n gmadet python=3 numpy scipy matplotlib astropy pandas shapely requests h5py scikit-image


Activate the environment:
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console
 
    $ conda activate gmadet 


Install other libraries
^^^^^^^^^^^^^^^^^^^^^^^

Once you have activated the environment, install the packages that are not available with conda using pip:

.. code-block:: console
 
    $ python3 -m pip install lacosmic hjson voevent-parse xmltodict astroML regions photutils keras keras-vis tensorflow cython regions  opencv-python-headless
    $ python3 -m pip install --pre astroquery

Install C dependencies:
^^^^^^^^^^^^^^^^^^^^^^^

* SExtractor: https://github.com/astromatic/sextractor
* SWarp: https://github.com/astromatic/swarp
* PSFEx: https://github.com/astromatic/psfex
* SCAMP: https://github.com/astromatic/scamp
* hotpants: https://github.com/acbecker/hotpants


.. _Github repo: https://github.com/dcorre/gmadet
.. _tarball: https://github.com/dcorre/gmadet/tarball/master



