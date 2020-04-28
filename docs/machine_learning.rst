=====================================
Machine learning to filter candidates
=====================================


Objective
---------

The image substraction results in many false candidates due to bad astrometric calibration, bad kernel estimation, bad photometric calibration, etc...


These artefacts are usually easy to discard by eye, but it is time consuming which is bad in an automatic process. So one soluton is to make use of machine learning to filter these false candidates, using a CNN algorithm. 

To do that, we run the pipeline on real images of a given telescope and classify by eye the true and false candidates. It means that we need to identify real transients or variable object in the images, which can be time consuming. We can also simulate point like objects in these astronomical images and run the pipeline on these images knowng exactly where these simulated transients are. To do that we estimate the PSF on each image, and use it to simulate some point like sources in these images. Then the pipeline is ran and the simulated transients are automatically classified as true candidates, the user has to examine by eye the rest of the candidates to determine whether they are true or false candidates. Then this set of true and false candidates are provided to the CNN algorithm to learn to classify candidates. Once the CNN training model is satisfactory, it can be used in the pipeline process to filter the candidates resulting of image substraction.


Process
-------

Run Docker if needed
^^^^^^^^^^^^^^^^^^^^

If you are using a Docker, run it first:

.. code-block:: console

    sudo docker run -v $(pwd)/gmadet/:/home/newuser/gmadet/ --rm -it gmadet:1.0

If the command $(pwd) does not work on your computer, just write the whole to the gmadet/gmadet/ repository on your computer.


Define test images 
^^^^^^^^^^^^^^^^^^

The first thing to do is to put some images that we will be used for simulating point like sources. Create a folder whose name corresponds to your telescope in /gmadet/gmadet/cnn/data/ . 


Estimate PSF of the test images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type: 

.. code-block:: console

    python psf.py --path data/yourtelescopeID/ --telescope yourtelescopeID

where 'yourtelescopeID' is the name of the folder where you put your images. The argument for --telescope is the telescope ID, corresponding to the folder name in gamdet/gmadet/config/ .

As a result, in 'data/yourtelescopename/' a '.psf.fits' file is created for each image. This file corresponds to the PSF estimated by PSFex, using the configuration parameters in gmadet/gmadet/config/ .


Simulate point like sources
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type:

.. code-block:: console

    python sim.py --datapath data/yourtelescopeID/ --telescope yourtelescopeID

The new fits images are created in gmadet/gmadet/cnn/sim/yourtelescopeID/images

A file named 'simulated_objects.list' will be also created, containing the positions of these new point like sources in the images.


Run gmadet
^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/ folder and type:

.. code-block:: console

    python run_gmadet.py --filename cnn/data/sim/yourtelescopeID/images/ --FWHM psfex --threshold 4 --radius_crossmatch 2.5 --telescope yourtelescopeID --doAstrometry scamp --doSub ps1



Create sub-images centered on each candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type:

.. code-block:: console

    python makesubimage.py --path cnn/data/sim/yourtelescopeID/images/ --telescope yourtelescopeID --training

The candidates are created in gmadet/gmadet/cnn/sim/yourtelescopeID/candidates/ . Two others folder are created true/ and false/, they will be used to contain what we classify as true and false candidates.
The simulated candidates are put in the true folder.

Classify true and false candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The idea is to create 2 folders, one for the true candidates and one for the false candidates. You can classify them by eye, perform a crossmatch with variable stars catalogs, etc...
The main thing is to put what you consider true and false candidates in the respective folders.


Run the CNN algorithm
^^^^^^^^^^^^^^^^^^^^^

Once you have classified your candidates, the next step is to trained the CNN algortihm to classify candidates. Before starting the training, we create a .npz datacube containing the candidates. Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type:

.. code-block:: console

    python convert.py --path data/sim/yourtelescopeID/candidates/ --telescope yourtelescopeID --cubename test

Then you can start the training:

.. code-block:: console

    python train.py --telescope yourtelescopeID --cubename test --modelname test


Use this model in gmadet
^^^^^^^^^^^^^^^^^^^^^^^^

