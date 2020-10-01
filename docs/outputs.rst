=======
Outputs
=======

(Need to be more detailed)

Short explanation of the different output files. The convention is that the original name of the image is kept and some keywords are appended to it or its extension.

* ``*_background``: background map used during the process.

* ``*_segmentation``: segmentation map provided by SExtractor.

* ``*_psf.fits``: PSF estimated in the image.

* ``*_SourcesDetnoFilter.cat``: raw output file of SExtractor.

* ``*_SourcesDet.cat``: reformated raw output file of SExtractor adding flag for bad sources.

* ``*.magwcs``: same as above but with RA, DEC position and a different format.

* ``*.magwcs2``: same as above but only with RA and DEC. Useful to load it in ds9 to visualise the sources.

* ``*.oc``: files listing all the sources in the original image not referenced in the catalogs.

* ``*.oc_RADEC``: same as above but containing the RA and DEC only. Useful to load them in ds9 for visualising.

* ``*_sub.oc`` and ``_sub.oc_RADEC``: same as above but for the substracted image.

* ``*.alldetections``: file containing the list of all sources in original and substracted images. Their photometry has been calibrated.

* ``*_transient_candidates.dat``: file containing the final list of possible transients. Their photometry has been calibrated. (Undergoing tests at the moment...)

* ``*_ZP_?_nosigmaClipping.png``: plot showing zeropoint without sigma clipping.

* ``*_ZP_?.png``: plot showing zeropoint with sigma clipping.

* ``*_ZP_?.dat``: file containing the result of the crossmatch between our sources and a catalog. Basically all the magnitudes in different bands and various flags.


* ``*_hotpants.sh``: commands used to run hotpants in a terminal. Useful when trying to tune hotpants parameters without re-running the whole process.

In case substraction was asked, a ``substraction/`` folder is created and additional files are present:


* ``*_mosaic``: original image recomposed when substraction was performed on each inidvidual PS1 images.

* ``*_mosaic_ps1``: mosaic of all the PS1 images matching the original image field of view.

* ``*_mosaic_sub``: substracted images of these mosiacs.

* ``*_mosaic_sub_mask``: mask applied during the substraction process.

* ``*_reg``: just means that images were registered.

* ``_0``, ``_1``, etc : index for each subimage that were created to perform a substraction on each individual PS1 image.



 





