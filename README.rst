======
gmadet
======


.. image:: https://img.shields.io/pypi/v/gmadet.svg
        :target: https://pypi.python.org/pypi/gmadet

.. image:: https://img.shields.io/travis/dcorre/gmadet.svg
        :target: https://travis-ci.org/dcorre/gmadet

.. image:: https://readthedocs.org/projects/gmadet/badge/?version=latest
        :target: https://gmadet.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Tools to help identification of transients in astronomical images. 


* Free software: MIT license
* Documentation: https://gmadet.readthedocs.io.


Features
--------

* Configurable for any telescopes through individual config files.
* Astrometric solution computed using SCAMP.
* Sources detection using Sextractor or pyRAF.
* Crossmatch of detected sources with Pan-STARRS DR1, GAIA DR2, GSC and USNO-B1 catalogs using the Xmatch algorithm.
* Crossmatch of detected sources with solar system moving objects using SkyBot.  
* Image substraction using Pan-STARRS stack images as reference:   
   * Image alignment performed using SWARP 
   * PSF estimation using PSFex
   * Substraction performed using hotpants
* Images can be stacked using SWARP before the analysis.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
