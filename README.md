# gmadet: GRANDMA detection pipeline

* Free software: MIT license
* Documentation: https://gmadet.readthedocs.io.

Development status
--------------------

[![Build Status](https://travis-ci.com/dcorre/gmadet.svg?branch=master)](https://travis-ci.com/dcorre/gmadet)
[![codecov](https://codecov.io/gh/dcorre/gmadet/branch/master/graphs/badge.svg)](https://codecov.io/gh/dcorre/gmadet/branch/master)
[![Documentation Status](https://readthedocs.org/projects/gmadet/badge/?version=latest)](https://gmadet.readthedocs.io/en/latest/?badge=latest)


Tools to help identification of transients in astronomical images. 

Aim
---

* Implementing a common detection pipeline for all the telescopes in the network.
* Each relevant telescope characteristics into account through individual configuration files for each external software being used (Astromatic softwares, hotpants).
* OS-independent, so that it can work on Windows, Mac OS and Linux machines.
* Written in Python3, with some external C dependencies.


Features
--------

* Usage of Docker to allow its deployment on different OS. Container is an Ubuntu 18.04 image. Successfully tested on Windows, Linux and Mac OS.
* Cosmics can be removed using the python library [lacosmic](https://github.com/larrybradley/lacosmic).
* Astrometric solution based on Gaia DR1 catalog computed using [SCAMP](https://github.com/astromatic/scamp).
* Image alignment, stacking, resampling using [SWarp](https://github.com/astromatic/swarp).
* Individual (part of) image PSF estimated with [PSFEx](https://github.com/astromatic/psfex).
* Source extractions using [SExtractor](https://github.com/astromatic/sextractor).
* Image substraction with Pan-STARRS reference images using [hotpants](https://github.com/acbecker/hotpants). PS1 images matching your image FoV are automatically downloaded.
* Crossmatch of detected sources with Pan-STARRS DR1, GAIA DR2, GSC and USNO-B1 catalogs using the Xmatch algorithm, through [astroquery](https://astroquery.readthedocs.io/en/latest/xmatch/xmatch.html).
* Crossmatch of detected sources with solar system moving objects using [SkyBot](https://astroquery.readthedocs.io/en/latest/imcce/imcce.html).
* Detected sources not referenced in any catalog are then classified using a CNN network implemented with [Keras](https://keras.io/, in order to keep only point-like sources).
* Photometric calibration (in AB) using Pan-STARRS DR1, SDSS DR12 or GAIA DR2 based on the zeropoint computed over the whole image.
* Automatic report of candidates to the database used for O3 using a HTTP POST method.


Installation
------------

See documentation: https://gmadet.readthedocs.io.


Credits
-------

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage) project template.
