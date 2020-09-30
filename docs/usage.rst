=====
Usage
=====

Testing that Docker image is working
------------------------------------

Run Docker
^^^^^^^^^^^^^^

Run the Docker image:

.. code-block:: console

   $ docker run -v /your_path_to_gmadet/:/home/newuser/gmadet/ -v /path_to_your_data/:/home/newuser/data/ --rm -it dcorre/gmadet

This means that you run interactively in a bash terminal the Docker image named dcorre/gmadet.  
The -v option means that you mount a volume in the Docker pointing to a directory on your computer. This allows to exchange data between the Docker and your machine. The first volume is pointing to the gmadet directory on your machine (the directory where the setup.py is). The second volume is pointing to the directory containing your images on your machine. For both cases, you need to edit the path before the `:`.


Install gmadet inside the Docker image.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   $ python3.7 setup.py develop --user

You will need to do it each time you run the Docker image, as everything is lost when you exit it (except the files created/edited in the mounted volumes).


Run gmadet on a test image.
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test if gmadet is running normally:

.. code-block:: console

    gmadet-run --path_data gmadet/data_test/ATLAS18qqn-S001-R001-C001-SDSS_g.fits --FWHM psfex --telescope IRIS --doAstrometry scamp --doSub ps1

It can take some times as it will download some Pan-STARRS archive image to perform the substraction. If it ran well you will see the last line starting with "Cleaning up output files for ...".


Description of the process
--------------------------

Functionalities
---------------
