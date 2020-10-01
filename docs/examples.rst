========
Examples
========

Run the Doker image
-------------------

If you are using the Docker image you need to run it first.

.. code-block:: console

   docker run -v /your_path_to_gmadet/:/home/newuser/gmadet/ -v /path_to_your_data/:/home/newuser/data/ --rm -it dcorre/gmadet

Replace:


* ``/your_path_to_gmadet/`` with the path on your machine pointing to the gamdet directory containing the ``setup.py``.   
* ``/path_to_your_data/`` with the path on your machine pointing to the data you want to analyse.

Then you need to install ``gmadet`` inside the docker image:

.. code-block:: console

   python3.7 setup.py develop --user && cd ..

This must be done each time you run the Docker image.


Astrometric calibration
-----------------------

.. code-block:: console

   gmadet-astrometry --path_data data/ --telescope your-telescope-alias

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-astrometry -h`` to see available ones.

This will perform the astrometric calibration on all the images contained in ``/data`` which is poiting to your data path defined above. If you want want to do it on a specific image or folder, specify it. The results will be stored in the folder ``gmadet_astrometry/``.

The ``--path_data`` can be a folder or a single image.

Type ``gmadet-astrometry -h`` to see the other arguments you might want to change.


Compute PSF
-----------

.. code-block:: console

   gmadet-psf  --path_data data/ --telescope your-telescope-alias

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-astrometry -h`` to see avai
lable ones.

This will compute the PSF using PSFEx for all images in ``/data``. The results will be stored in the folder ``gmadet_psf/``.

Type ``gmadet-psf -h`` to see the other arguments you might want to change.


Stacking
--------

.. code-block:: console

   gmadet-stacking --path_data data/ --radius RADIUS --deltaT DELTAT

Replace:

* ``RADIUS``: with a float with units arcminutes. All images whose center is below this value will be considered as the same field.   
* ``DELTAT``: with a float with units hour. All images of the same field taken with the same telescope and instrument and filter band taken within this duration will be stacked.

The results will be stored in the folder ``gmadet_stacking/``.


Substracting background
-----------------------

.. code-block:: console

   gmadet-subBkg --path_data data/ 

This will substract the background of all images in ``data/`` using the same method as SExtractor by default. Type ``gmadet-subBkg -h`` to see the other arguments you might want to change. The results are stored in ``gmadet_subBkg/``.


Remove cosmics
--------------

.. code-block:: console

   gmadet-cosmics --path_data data/

This will remove cosmic rays using the L.A. Cosmic algorithm. Results are stored in ``gmadet_remove_cosmics/``.

Type ``gmadet-cosmics -h`` to see the other arguments you might want to change.

Following the documentation, 4 iterations should be the maximum, if sources are still removed after you are likely removing pixels from saturated stars for instance.


Run gmadet without image substraction
-------------------------------------

.. code-block:: console

   gmadet-run --path_data data/ --FWHM psfex --telescope your-telescope-alias --doAstrometry scamp --radius_crossmatch 3 --threshold 4

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-run -h`` to see available ones.

For all images in ``data/`` this will perform:

* Astrometric calibration with SCAMP using GAIA DR2 by default.   
* Find sources using SExtractor using a threshold of 4.  
* PSFEx is sued to estimate the PSF FWHM.    
* Crossmatch all sources with catalogs (GAIA DR2, PS1 DR1, GSC, USNO-B1) within 3 pixels. Xmatch is used to do the crossmatch with online queries.   
* Crossmatch solar moving objects using SkyBoT.   

Type ``gmadet-run -h`` to see the other arguments you might want to change. You can add backgroung subtraction, removal of cosmics for instance.

Results are stored in ``gmadet_results/``.

Run gmadet with image substraction using PS1 image reference
------------------------------------------------------------

.. code-block:: console

   gmadet-run --path_data data/ --FWHM psfex --telescope your-telescope-alias --doAstrometry scamp --radius_crossmatch 3 --threshold 4 --doSub ps1 --ps1_method individual

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-run -h`` to see available ones.

For all images in ``data/`` this will perform:

* Astrometric calibration with SCAMP using GAIA DR2 by default.   
* PSFEx is sued to estimate the PSF FWHM.    
* If not already present in ``gmadet/ps1Dir/``, download PS1 archive stack images matching your image field of view. Then rescale the images to a linear scale and store them in ``gmadet/ps1RescaledDir/``.   
* Perform an image substraction using hotpants. The ``--ps1_method individual`` means that the substraction will be performed using each PS1 images separately. All subimages are combined in a substracted mosaic image at the end of the process.   
* Find sources using SExtractor on the substracted mosaic image using a threshold of 4.   
* Crossmatch all sources with catalogs (GAIA DR2, PS1 DR1, GSC, USNO-B1) within 3 pixels. Xmatch is used to do the crossmatch with online queries.   
* Crossmatch solar moving objects using SkyBoT.   

Type ``gmadet-run -h`` to see the other arguments you might want to change. You can add backgroung sub
traction, removal of cosmics for instance.

Results are stored in ``gmadet_results/``. Result of substraction in ``gmadet_results/substraction/``.


**IMPORTANT**: PS1 survey is limited to -30 degrees in declination, so can only be used above.

