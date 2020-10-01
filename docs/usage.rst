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
The -v option means that you mount a volume in the Docker pointing to a directory on your computer. This allows to exchange data between the Docker and your machine. The first volume is pointing to the gmadet directory on your machine (the directory where the setup.py is). The second volume is pointing to the directory containing your images on your machine. For both cases, you need to edit the path before the ``:``.


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

Telescopes
^^^^^^^^^^

The code is configurable for each telescope through configuration files for each software (SExtractor, SWarp, PSFEx, SCAMP, hotpants). These configuration files are in ``gmadet/config/``. You need to edit these files to tune the code for your telescope.

To add one telescope, simply add a folder with a telescope alias in a similar way. You can copy an existing one and then edit the files and change the name of the folder.

To specify which telescope must be used, you will provide the alias of the telescope to the executables through the ``--telecope`` argument. The alias is simply the name of the folder in ``gmadet/config/``.

Executables
^^^^^^^^^^^

All the scripts can be run with executable beginning with "gamdet-". Then you can provide the options through arguments from the command line, for instance type ``gmadet-astrometry -h`` to know the expected arguments for the executable performing astrometry calibration. Below the list of executables:

* **gmadet-astrometry**: perform astrometry usig SCAMP. 

* **gmadet-psf**: Compute the PSF using PSFEx.

* **gmadet-stacking**: Stack images using SWarp.

* **gmadet-subBkg**: Substract background using the some of the routines of `photutils <https://photutils.readthedocs.io/en/stable/background.html>`.

* **gmadet-cosmics**: Remove cosmics rays using `L.A. Cosmic algorithm <https://lacosmic.readthedocs.io/en/latest/>`.



Functionalities
---------------
