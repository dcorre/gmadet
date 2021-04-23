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

    git clone git://github.com/dcorre/gmadet

Note that you will need to update the project regularly to check for updates. Ideally each time you want to use it, type the following command to search for updates:

.. code-block:: console

    git pull origin master


* Or download:

  Simply download this `tarball`_. Or through the console:

.. code-block:: console

    curl  -OL https://github.com/dcorre/gmadet/tarball/master

Cloning the project allows to retrieve easily future code updates. If you downloaded the project, you will need to download it again to retrieve future updates.


Installation with Docker
------------------------

The usage of a Docker allows to build an OS environment on your machine and thus avoid compatibility problems when running the code under Linux, Mac or Windows. If you have not Docker installed on your machine install it first.

* Install the Docker desktop for your OS: https://docs.docker.com/get-docker/

* To run Docker without appending sudo, type:

.. code-block:: console

   sudo groupadd docker
   sudo usermod -aG docker $USER

Log out and log back in so that your group membership is re-evaluated. For more information see https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user.

You can test that Docker is installed correctly and can be run without sudo:

.. code-block:: console

   docker run hello-world


Download the gmadet Docker image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To retrieve the Docker image:

.. code-block:: console

   docker pull dcorre/gmadet

Check that it appears in the list of images:

.. code-block:: console

   docker images


Installation without Docker
---------------------------

I advise to create a virtual environment to avoid messing with different python libraries version that could be already installed on your computer and required for other projects.

Install conda: https://docs.conda.io/en/latest/miniconda.html

You can also install everything with pip if you prefer not to use conda.

Python 3 environment:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    conda create -n gmadet python=3 numpy scipy matplotlib astropy pandas shapely requests h5py scikit-image


Activate the environment:
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    conda activate gmadet


Install other libraries
^^^^^^^^^^^^^^^^^^^^^^^

Once you have activated the environment, install the packages that are not available with conda using pip:

.. code-block:: console

    python3 -m pip install lacosmic hjson voevent-parse xmltodict astroML regions photutils keras keras-vis tensorflow cython regions opencv-python-headless astroquery

Install C dependencies:
^^^^^^^^^^^^^^^^^^^^^^^

* SExtractor: https://github.com/astromatic/sextractor
* SWarp: https://github.com/astromatic/swarp
* PSFEx: https://github.com/astromatic/psfex
* SCAMP: https://github.com/astromatic/scamp
* hotpants: https://github.com/acbecker/hotpants


.. _Github repo: https://github.com/dcorre/gmadet
.. _tarball: https://github.com/dcorre/gmadet/tarball/master


Testing that it is working
--------------------------

Run Docker
^^^^^^^^^^^^^^

Run the Docker image in the background:

.. code-block:: console

   docker run --name gmad -dit -v /your_path_to_gmadet/:/home/newuser/gmadet/ -v /path_to_your_data/:/home/newuser/data/ dcorre/gmadet

| This means that you run the docker image `dcorre/gmadet`, give the name `gmad` to the created container.
| `-d` runs the container in backgound.   
| `-i` gives the possibility to enter in the container to run commands interactively in a bash terminal.
| `-t` allocates a pseudo-TTY. 
| The -v option means that you mount a volume in the Docker pointing to a directory on your computer. This allows to exchange data between the Docker and your machine.
| The first volume is pointing to the gmadet directory on your machine (the directory where the setup.py is). The second volume is pointing to the directory containing your images on your machine. For both cases, you need to edit the path before the ``:``.

Once you have executed this command, you can run any command in the container by typing:

.. code-block:: console

   docker exec gmad ls
   docker exec gmad pwd
   
to make a `ls` or a `pwd` in the container named `gmad`, or any other bash commands.

**In the following, if you are using a Docker image just prepend the command `docker exec gmad` to run the given commands within the container instead of your machine.** 

The container is alive as long as you do not shut down your machine. It is important to know that you can not give the same name to two containers. So if for some reasons you need to remove the current container to start a new one, type:

.. code-block:: console

   docker rm gmad

You can list the containers, active or not, on your machine with:

.. code-block:: console

   docker ps -a



Install gmadet
^^^^^^^^^^^^^^
-------------------------
Inside the Docker image
-------------------------

.. code-block:: console

   docker exec gmad python3.8 setup.py develop --user

It is valid as long as the container is alive. You will need to do it each time you start a new container.  the Docker image, as everything is lost when you exit it (except the files created/edited in the mounted volumes).

--------------------
Without Docker image
--------------------

.. code-block:: console

   python3 setup.py develop



Run gmadet on a test image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test if gmadet is running normally:

.. code-block:: console

   gmadet-run gmadet/data_test/test_image.fits --telescope IRIS --sub ps1 --result gmadet/data_test/gmadet_results

It can take some times as it will download some Pan-STARRS archive image to perform the substraction. If it ran well you will see the last line starting with "Cleaning up output files for ...".
A folder gmadet_results/test_image has been created in gmadet/data_test/ with a bunch of files that will be described later on.
