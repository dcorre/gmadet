=====================================
Machine learning to filter candidates
=====================================


Objective
---------

The image substraction results in many false candidates due to bad astrometric calibration, bad kernel estimation, bad photometric calibration, etc...


These artefacts are usually easy to discard by eye, but it is time consuming and obviously not possible in an automatic process. So one soluton is to make use of machine learning to filter these false candidates, using for instance a CNN algorithm.

To do that, we can run the pipeline on real images of a given telescope and classify by eye the true and false candidates. It means that we need to identify real transients or variable object in the images, which can be time consuming.

We can also simulate point like objects in these astronomical images and run the pipeline on these images knowing exactly where these simulated transients are. To do that we estimate the PSF on each image, and use it to simulate some point like sources in these images. Then the pipeline is ran and the simulated transients are automatically classified as true candidates, the user can either examine by eye the rest of the candidates to determine whether they are true or false candidates or consider all of them as false candidates (not the best but if the number of simulated sources is large enough it can be a good comprise). Then this set of true and false candidates are provided to the CNN algorithm to learn to classify candidates. Once the CNN training model is satisfactory, it can be used in the pipeline process to filter the candidates resulting of image substraction.


Process
-------

Run the Doker image
^^^^^^^^^^^^^^^^^^^

If you are using the Docker image you need to run it first.

.. code-block:: console

   docker run -v /your_path_to_gmadet/:/home/newuser/gmadet/ -v /path_to_your_data/:/home/newuser/data/ --rm -it dcorre/gmadet

Replace:
* ``/your_path_to_gmadet/`` with the path on your machine pointing to the gamdet directory containing the ``setup.py``.
* ``/path_to_your_data/`` with the path on your machine pointing to the data you want to use for the CNN training.

Then you need to install ``gmadet`` inside the docker image:

.. code-block:: console

   python3.7 setup.py develop --user && cd ..

This must be done each time you run the Docker image.



Simulate point-like sources in your images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type: 

.. code-block:: console

   gmadet-sim --path_data data/ --telescope your-telescope-alias --Ntrans 100 --magrange 14 23 --ZP 30

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-astrometry -h`` to see avai
lable ones.

This will insert 100 sources on each of your image, taking the PSF of each image (or part of image), whose magnitude in the range 14-23 with a zeropoint of 30.

Results are stored in ``gmadet_sim/simulation/``. A file named 'simulated_objects.list' will be also created, containing the positions of these new point like sources in the images.


Type ``gmadet-sim -h`` to see the other arguments you might want to change.


Run gmadet on these simulated images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We run it with the image substraction but you can do it without.

.. code-block:: console

   gmadet-run --path_data data/gmadet_sim/simulation/ --FWHM psfex --telescope your-telescope-alias --doAstrometry scamp --radius_crossmatch 3 --threshold 4 --doSub ps1 --ps1_method individual

The results are stored in ``data/gmadet_sim/simulation/gmadet_results``.


Create image cutouts of all candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    gmadet-cutouts --path_data data/gmadet_sim/simulation/ --training 

The ``--training`` argument specifies that it is for the training on simulated images and create a ``true`` and ``false`` folder in ``candidates_training``. They will be used for the CNN training as what we consider true and false candidates. The simulated candidates are automatically put in the ``true`` folder.

You can either classify the other ones by eye, or put all of them in the ``false`` folder. Some true sources will be classified as false but if the number of simulated sources is large enough, this might be a comprise.


You can plot some histograms to check the distribution of magnitudes for the different bands and fraction of the simulated objects that are actually detected by writing:

.. code-block:: console

    gmadet-checksim --path_data data/gmadet_sim/simulation/  --radius 2

It will create a folder ``CheckSim/`` with some plots. It will also create a file ``crossmatch.dat`` containing the crossmatch of the sources detected by gmadet and the positions of the simulated sources. This is useful to make some tests of how the code behaves with known simulated transients.


Classify true and false candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The idea is to create 2 folders, one for the true candidates and one for the false candidates. You can classify them by eye, perform a crossmatch with variable stars catalogs, etc...
The main thing is to put what you consider true and false candidates in the respective folders.


Run the CNN algorithm
^^^^^^^^^^^^^^^^^^^^^

Once you have classified your candidates, the next step is to trained the CNN algorithm to classify candidates. Before starting the training, we need to create a .npz datacube containing the candidates in the right format.

.. code-block:: console

    gmadet-cnn_convert --path_datacube PATH_DATACUBE --cubename CUBENAME --path_cutouts PATH_CUTOUTS

Replace:

* ``PATH_DATACUBE`` with the pah where you want to store your datacube. 
* ``CUBENAME`` with the name of the datacube that will be created.
* ``PATH_CUTOUTS`` with the path to the folder containing the ``true`` and ``false`` folders.

Then you can start the training:

.. code-block:: console

    gmadet-cnn_train --path_cubename PATH_CUBENAME --path_model PATH_MODEL --modelname MODELNAME

Replace:

* ``PATH_CUBENAME`` with the path containing the datacube, including the filename and .npz extension.
* ``PATH_MODEL`` with the path where you want to store the trained model.
* ``MODELNAME`` with the name of the model that will be created.


You can visualise the results with some plots that will help to assess the probability threshold to apply.

.. code-block:: console

    gmadet-cnn_checkinfer --path_plots PATH_PLOTS --path_crossmatch PATH_CROSSMATCH --path_infer PATH_INFER 

Replace:

* ``PATH_PLOTS`` with the path where you want to store the plots.
* ``PATH_CROSSMATCH`` with the path where the crossmatch.dat is stored.
* ``PATH_INFER`` with the path where the infer_results.dat is stored.


Type ``gmadet-cnn_checkinfer -h`` to see the other optional arguments.

Apply a trained model on candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It assumes that you already ran gmadet on a set of images, and created the candidates cutouts without using the ``--training`` argument. So you will have a ``candidates`` folder containing the cutouts that need to be classify by the CNN algorithm.

.. code-block:: console

    gmadet-cnn_infer --path_cutouts PATH_CUTOUTS --path_model PATH_MODEL

Replace:

* ``PATH_CUTOUTS`` with the path containing the candidates cutouts.
* ``PATH_MODEL`` with the path to the trained CNN model, including its filnemame and .h5 extension.

It will results a file ``infer_results.dat`` in the ``PATH_CUTOUTS`` containing the probability that a source is a false (column: label0) or true (column: label1). You can then aplly a threshold on these probability to keep only some candidates.

To assess the threshold you can run the ``gmadet-cnn_infer`` on the same images you used for the training.

Then you can visualise the results with some plots that will help to assess the probability threshold to apply.

.. code-block:: console

    gmadet-cnn_checkinfer --path_plots PATH_PLOTS --path_crossmatch PATH_CROSSMATCH --path_infer PATH_INFER 

Replace:

* ``PATH_PLOTS`` with the path where you want to store the plots.
* ``PATH_CROSSMATCH`` with the path where the ``crossmatch.dat`` is stored.
* ``PATH_INFER`` with the path where the ``infer_results.dat`` is stored.


Type ``gmadet-cnn_checkinfer -h`` to see the other optional arguments.

