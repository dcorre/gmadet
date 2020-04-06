=====
Usage
=====

If you are using a Docker, run it first:

.. code-block:: console

    sudo docker run -v $(pwd)/gmadet/:/home/newuser/gmadet/ --rm -it gmadet:1.0

To test if gmadet is running normally, go to gmadet directory where gmadet.py is and type::

    python gmadet.py --filename test/ATLAS18qqn-S001-R001-C001-SDSS_g.fits --FWHM psfex --threshold 4 --radius_crossmatch 2.5 --telescope IRIS --doAstrometry scamp --doSub ps1
