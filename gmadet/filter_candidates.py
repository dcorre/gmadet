#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import join
from gmadet.cnn.infer import infer
from gmadet.utils import (make_sub_image, mkdir_p)


def filter_candidates(sources,
                      FWHM_ratio_lower=0.5,
                      FWHM_ratio_upper=5.0,
                      CNN_model=None,
                      CNN_thres=0.8,
                      makecutout=True,
                      size=32,
                      fmt='png',
                      outLevel=1):
    """Filter transient candidates"""
    print ('Filter candidates')

    # Take first candidate to extract the path where to store the file
    # No need to chack if substraction was performed, 
    # as if it did only the ones from substracted files are with 'Match' == Y
    path, fname_ext = os.path.split(sources['filenames'][0])
    if path:
        path = path + "/"
    # Get rid of the extension to keep only the name
    fname2, extension = os.path.splitext(fname_ext)
    # Get rid of the _reg pattern
    fname2 = fname2.split('_ref')[0]

    # First get the sources not crossmatching with sources in catalogs
    mask_cat = sources['Match'] == 'N'

    # Remove candidates on the edges
    mask_edge = sources["edge"] == 'N'

    # Remove sources with FWHM ratio outside the desired range
    FWHM_ratio = sources["FWHM"] / sources["FWHMPSF"]
    mask_FWHM = (FWHM_ratio >= FWHM_ratio_lower) & \
                (FWHM_ratio <= FWHM_ratio_upper)

    mask_tot = mask_cat & mask_edge & mask_FWHM

    # Use a trained CNN model to filter candidates.
    if CNN_model is not None:
        # Create fits cutouts to be be given to the CNN model
        path_CNN_cutout = path + fname2 + '_CNN_cutouts/'
        mkdir_p(path_CNN_cutout)
        for cand in candidates:
            coords = [cand['_RAJ2000'], cand['_DEJ2000']]
            outname = path_CNN_cutout + 'candidate_%d.%s' % (cand['cand_ID'],
                                                             fmt)
            info_dict = {}
            info_dict['RA'] = cand['_RAJ2000'] 
            info_dict['DEC'] = cand['_DEJ2000']
            info_dict['XPOS'] = cand['Xpos']
            info_dict['YPOS'] = cand['Ypos']
            info_dict['FILE'] = cand['filenames']
            info_dict['CANDID'] = cand['cand_ID']
            info_dict['MAG'] = cand['mag_calib']
            info_dict['MAGERR'] = cand['mag_calib_err']
            info_dict['FWHM'] = cand['FWHM']
            info_dict['FWHMPSF'] = cand['FWHMPSF']
            make_sub_image(
                cand['filenames'],
                coords,
                coords_type="world",
                output_name=outname,
                size=[size, size],
                fmt="fits",
                addheader=True,
                title=None,
                info_dict=info_dict
            )
        # Run CNN model to associate a probability to each cutout
        # The size of the cutout should be the same as the ones used
        # for the CNN training
        infer(path_CNN_cutouts, CNN_model, 0.0)

        # Add the probability to the canidates table.
        infer_table = ascii.read(path_CNN_cutout + 'infer_results.dat')
        candidates = join(candidates,
                         infer_table['cand_ID', 'label0', 'label1'],
                         join_type='left')

        # keep only transients that are above the threshold
        mask_CNN = candidates['label1'] >= CNN_thres
        mask_tot = mask_tot & mask_CNN

    # Write output file.
    candidates = sources[mask_tot]
    # Create ID to start from 1.
    candidates['cand_ID'] = np.arange(len(candidates)) + 1
    # Rename colums
    if 'label0' in candidates.colnames:
        candidates.rename_column('label0', 'P_False')
        candidates.rename_column('label1', 'P_True')
    candidates.write(
        path + fname2 + '_candidates.dat',
        format='ascii.commented_header',
        overwrite=True
    )

    # Update 
    print ('Make cutouts')
    #  Extract small image centered on candidates passing the filters
    if makecutout:
        path_cutout = path + fname2 + '_cutouts/'
        mkdir_p(path_cutout)
        for cand in candidates:
            coords = [cand['_RAJ2000'], cand['_DEJ2000']]
            outname = path_cutout + 'candidate_%d.%s' % (cand['cand_ID'],
                                                         fmt)
            if fmt == 'png':
                _coords = SkyCoord(cand['_RAJ2000'], cand['_DEJ2000'],
                                  unit=(u.degree,u.degree),
                                  frame='icrs')
                coords_sexa = _coords.to_string(style='hmsdms')
                title = 'RA Dec: %s \n Mag: %.2f +/- %.2f \nFWHM_ratio: %.2f' % \
                        (coords_sexa,
                         cand['mag_calib'],
                         cand['mag_calib_err'],
                         cand['FWHM']/cand['FWHMPSF'])
            make_sub_image(
                cand['filenames'],
                coords,
                coords_type="world",
                output_name=outname,
                size=[size, size],
                fmt=fmt,
                addheader=False,
                title=title
            )


