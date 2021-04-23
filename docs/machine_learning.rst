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

Launch the Doker image
^^^^^^^^^^^^^^^^^^^



If you are using the Docker image, remember to launch once the container:

.. code-block:: console

   docker run --name gmad -dit -v /your_path_to_gmadet/:/home/newuser/gmadet/ -v /path_to_your_data/:/home/newuser/data/  dcorre/gmadet

Replace:


* ``/your_path_to_gmadet/`` with the path on your machine pointing to the gamdet directory containing the ``setup.py``.
* ``/path_to_your_data/`` with the path on your machine pointing to the data you want to analyse.

Then you need to install ``gmadet`` inside the docker image:

.. code-block:: console

   python3.8 setup.py develop --user

This must be done each time you launch the Docker image.

Then you only need to prepend `docker exec gmad` to the commands given below to execute them within the container instead of your machine.


Simulate point-like sources in your images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whithin a terminal, go to the gmadet/gmadet/cnn/ folder and type:

.. code-block:: console

   gmadet-sim data/  --results gmadet_sim --telescope your-telescope-alias --ntrans 100 --magrange 14 23 --zp 30

Replace:

* ``your-telescope-alias`` with the alias of your telescope, type ``gmadet-astrometry -h`` to see available ones.

This will insert 100 sources on each of your image, taking the PSF of each image (or part of image depending on your PSFex configuration file), whose magnitudes lie in the range 14-23 with a zeropoint of 30.

Results are stored in ``gmadet_sim/``. A file named `simulated_objects.list` is also created, containing the positions of these new point like sources in the images.


Type ``gmadet-sim -h`` to see the other arguments you might want to change.


Run gmadet on these simulated images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We run it with the image substraction but you can do it without.

.. code-block:: console

   gmadet-run gmadet_sim/ --results gmadet_sim_results --telescope your-telescope-alias --radius-crossmatch 3 --threshold 4 --sub ps1 --ps1-method individual

The results are stored in ``gmadet_sim_results``. Note that each file `simulated_objects.list` is copied to this new folder structure and the filenames defined in the 'filename' column are updated to the new location.


Create image cutouts of all candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    gmadet-cutouts gmadet_sim_results/ --training

The ``--training`` argument specifies that it is for the training on simulated images and create a ``true`` and ``false`` folders in ``candidates_training``. They will be used for the CNN training as what we consider true and false candidates. The simulated candidates are automatically put in the ``true`` folder.

You can either classify the other ones by eye, or put all of them in the ``false`` folder (use argument ``--false`` to do it automatically). Some true sources will be classified as false but if the number of simulated sources is large enough, this might be a comprise.


You can plot some histograms to check the distribution of magnitudes for the different bands and fraction of the simulated objects that are actually detected by writing:

.. code-block:: console

    gmadet-checksim gmadet_sim_results

It will create a folder ``CheckSim/`` with some plots. It will also create a file ``crossmatch.dat``, if not already created by ``gmadet-cutouts``, containing the crossmatch of the sources detected by gmadet and the positions of the simulated sources. This is useful to make some tests of how the code behaves with known simulated transients.


Classify true and false candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The idea is to create 2 folders, one for the true candidates and one for the false candidates. You can classify them by eye, perform a crossmatch with variable stars catalogs, etc...
The main thing is to put what you consider true and false candidates in the respective folders.


Run the CNN algorithm
^^^^^^^^^^^^^^^^^^^^^

Once you have classified your candidates, the next step is to trained the CNN algorithm to classify candidates. Before starting the training, we need to create a .npz datacube containing the candidates in the right format.

.. code-block:: console

    gmadet-cnn_convert --path PATH_DATACUBE --cube CUBENAME --cutouts PATH_CUTOUTS

Replace:

* ``PATH_DATACUBE`` with the pah where you want to store your datacube.
* ``CUBENAME`` with the name of the datacube that will be created.
* ``PATH_CUTOUTS`` with the path to the folder containing the ``true`` and ``false`` folders.

For the setup used in the previous examples, it will be

.. code-block:: console

    gmadet-cnn_convert --path gmadet_cnn --cube cube --cutouts gmadet_sim_results/candidates_training/

The cube will be in ``gmadet_cnn/datacube/cube.npz``

Then you can start the training:

.. code-block:: console

    gmadet-cnn_train --cube PATH_CUBENAME --model-path PATH_MODEL --model-name MODELNAME

Replace:

* ``PATH_CUBENAME`` with the path containing the datacube, including the filename and .npz extension.
* ``PATH_MODEL`` with the path where you want to store the trained model.
* ``MODELNAME`` with the name of the model that will be created.

Again, it will look like that

.. code-block:: console

    gmadet-cnn_train --cube gmadet_cnn/datacube/cube.npz --model-path gmadet_cnn --model-name model

The model will be in ``gmadet_cnn/CNN_training/model.h5``


Apply a trained model on candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It assumes that you already ran gmadet on a set of images, and created the candidates cutouts without using the ``--training`` argument. So you will have a ``candidates`` folder containing the cutouts that need to be classify by the CNN algorithm. 


.. code-block:: console

    gmadet-cnn_infer --cutouts PATH_CUTOUTS --model PATH_MODEL

Replace:

* ``PATH_CUTOUTS`` with the path containing the candidates cutouts.
* ``PATH_MODEL`` with the path to the trained CNN model, including its filnemame and .h5 extension.

For our example above, you can simply create a ``candidates`` folder in ``gmadet_sim_results/`` containing all the cutouts in ``gmadet_sim_results/candidates_training/`` true/ and false/ folders. This avoids re-running gmadet on the images, as we already did it for the training. Then we can apply the CNN trained model on the same cutouts we used for the training. If the training went well, it should classify all the simulated sources as real transients, apart from some of the faintest ones.

Again, for our eaxmaple it will look like that

.. code-block:: console

    gmadet-cnn_infer --cutouts gmadet_sim_results/candidates/ --model gmadet_cnn/CNN_training/model.h5

It will result a file ``infer_results.dat`` in the directory defined with ``--cutouts``, ``gmadet_sim_results/candidates/`` for our example, containing the probability that a source is a false (column: label0) or true (column: label1) transient.    
You can then apply a threshold on these probability to keep only some candidates. An idea would be to select the threshold according to the False Positive Rate, i.e. you select the probability corresponding to at most 1% (or whatever value suitable for you) of false positive in your trained sample.

To visualize how these probabilities evolve with some of the candidates parameters (magnitude, FWHM) of your sample, you can use ``gmadet-cnn_checkinfer``.

.. code-block:: console

    gmadet-cnn_checkinfer --plots PATH_PLOTS --crossmatch PATH_CROSSMATCH --infer PATH_INFER

Replace:

* ``PATH_PLOTS`` with the path where you want to store the plots.
* ``PATH_CROSSMATCH`` with the path where the ``crossmatch.dat`` is stored.
* ``PATH_INFER`` with the path where the ``infer_results.dat`` is stored.


Type ``gmadet-cnn_checkinfer -h`` to see the other optional arguments.

Again, for our example it will be

.. code-block:: console

    gmadet-cnn_checkinfer --plots gmadet_sim_results/ --crossmatch gmadet_sim_results/ --infer gmadet_sim_results/candidates/

It will results a folder ``CheckInfer`` in ``gmadet_sim_results/`` containing some plots illustrating the dependence of the probability that a candidate is a true transient (returned by the CNN algorithm) as a function of magnitude and FWHM ratio (so far, can include more check in the future). It also compares this evolution for the simulated soures with respect to the non-simulated sources. It is also useful to get an idea of the FWHM ratio range that can be applied to filter the candidates.

General notes
^^^^^^^^^^^^^

Ideally the training should be done on a few tens of images with taken in different observing conditions (elevation, seeing, moon phase, etc...) so that you can train a model that is representative enough of the images you can have, and thus not having to train a model for each sample of images you want to analyse.

Of course, if the computational time is not a constraint for you, it will be more accurate to perform a training on the images you want to analyse only, if you have a sufficient number of them.

Regarding the total number of transients required for an accurate training, I do not have a proper answer to that question. I would say the more the better. A few thousands to a few tens of thousands true transients should be enough, which is easy to achieve using the simulated sources. If you work with non-simulated sources, a few thousands of visually inspected objects might be enough, but again it is just a guess and the best thing is to try and see how the code behaves regarding to your scientific case.

Having a similar number of true and false transients in your training sample seems reasonable although I haven't tested the influence of having more false transients then true ones.
