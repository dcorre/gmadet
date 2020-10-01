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

All the scripts can be run with executable beginning with "gamdet-". Then you can provide the options through arguments from the command line, for instance type ``gmadet-astrometry -h`` to know the expected arguments for the executable performing astrometry calibration. Below the list of executables (use the ``-h`` options to see the expected arguments):

Main executables:

* **gmadet-astrometry**: perform astrometry usig SCAMP. 

* **gmadet-psf**: Compute the PSF using PSFEx.

* **gmadet-stacking**: Stack images using SWarp.

* **gmadet-subBkg**: Substract background using the some of the routines of `photutils`_.

* **gmadet-cosmics**: Remove cosmics rays using `L.A. Cosmic algorithm`_.

* **gmadet-run**: Main executable allowing to do all the taks described above, plus perform image substraction, crossmatch candidates with existing catalogs, perform photometric calibration. Can also apply a trained CNN model to classify transients. Results can be automatically reported to a database.

Executables related to the Convolutional Neural Network usage:

* **gmadet-sim**: Simulate point-like sources in images using the PSF of each image (or part of image) computing with PSFEx.

* **gmadet-cutouts**: Create image cutouts centered on the position of candidates.

* **gmadet-checksim**: Do some plots to visualise distribution of simulated sources for instance.

* **gmadet-cnn_convert**: Convert the candidates cutouts (classified as true/false events) into a single datacube with the format expected by the CNN algorithm.

* **gmadet-cnn_train**: Train the CNN algorithm with the created datacube using Keras.

* **gmadet-cnn_infer**: Apply a trained CNN model on a set of candidates cutouts to assign a probability of being a true or false event.

* **gmadet-cnn_checkinfer**: Do some plots to visualise the CNN training. Useful for chosing the probability threshold that will be used to classify a candidate as true or false event.

.. _photutils: https://photutils.readthedocs.io/en/stable/background.html
.. _L.A. Cosmic algorithm: https://lacosmic.readthedocs.io/en/latest/

Important information
---------------------

Image header
^^^^^^^^^^^^

As the headers can be heterogeneous, the code will only keep standard keywords and remove all the others. The kept keywords are:

* SIMPLE
* BITPIX
* NAXIS
* NAXIS1
* NAXIS2
* EXTEND
* DATE-OBS
* TELESCOP
* INSTRUME
* OBJECT
* EXPTIME
* FILTER
* GAIN
* SATURATE
* EQUINOX
* EPOCH
* RADESYS
* CTYPE1
* CTYPE2
* CUNIT1
* CUNIT2
* CRVAL1
* CRVAL2
* CRPIX1
* CRPIX2
* CD1_1
* CD1_2
* CD2_1
* CD2_2
* CDELT1
* CDELT2
* CROTA1
* CROTA2
* BSCALE
* BZERO
* BUNIT
* AIRMASS
* END


Check their definition `here`_

.. _here: https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html

To update the list of keyword, need to edit ``gmadet/sanitise.py``.


Astrometric calibration
^^^^^^^^^^^^^^^^^^^^^^^

The astrometric calibration is performed using SCAMP that requires an initial astrometric solution. Otherwise it can take a long time to run, and will most likely fails. The kept header keywords described above should be sufficient for SCAMP to converge rapidly.

Configuration files
^^^^^^^^^^^^^^^^^^^

You can find the description of the parameters used by the astromatic software and hotpants here:

* SExtractor:  https://sextractor.readthedocs.io/en/latest/ and http://www.astromatic.net/software/sextractor
* SWarp: http://www.astromatic.net/software/swarp
* PSFEx: https://psfex.readthedocs.io/en/latest/ and http://www.astromatic.net/software/psfex
* SCAMP: https://scamp.readthedocs.io/en/latest/ and http://www.astromatic.net/software/scamp
* hotpants: https://github.com/acbecker/hotpants

The default configuration can work reasonably well. 

It is important to note that the following parameter will be overwritten by the command line arguments:

**SExtractor:**

* ``DETECT_THRESH``: overwritten by the ``--threshold`` argument of ``gmadet-run``.
* ``FILTER_NAME``: overwritten by the ``--convFilter`` argument of ``gmadet-run`` and other executables.
* ``SEEING_FWHM``: overwritten by the ``--FWHM`` argument (either a float or the value returned by PSFEx)
* ``VERBOSE_TYPE``: overwritten by the ``--verbose`` argument of different executables.
* ``PARAMETERS_NAME``, ``CHECKIMAGE_TYPE``, ``CHECKIMAGE_NAME``, ``CATALOG_NAME``, ``PSF_NAME`` are also overwritten in the process.

**SCAMP:**

* ``FILTER_NAME``: overwritten by the ``--convFilter`` argument.
* ``VERBOSE_TYPE``: overwritten by the ``--verbose`` argument of different executables.
* ``PARAMETERS_NAME``, ``-ASTREF_BAND``, ``CHECKPLOT_DEV``, ``CHECKPLOT_NAME``, ``CHECKPLOT_TYPE`` are also overwritten in the process.

**PSFEx**:

* ``VERBOSE_TYPE``: overwritten by the ``--verbose`` argument of different executables.

For PSFEx, if you have a large field of view or non negligible spatial variation of the PSF, ``PSFVAR_NSNAP`` parameter controls how many PSF per axis you want to compute.
