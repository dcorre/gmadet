#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import join
from astropy.time import Time
from astropy.io import ascii, fits
from gmadet.cnn.infer import infer
from gmadet.utils import (
    make_sub_image,
    make_fits,
    make_figure,
    combine_cutouts,
    mkdir_p)
import multiprocessing as mp


def filter_candidates(sources,
                      FWHM_ratio_lower=0.5,
                      FWHM_ratio_upper=5.0,
                      CNN_model=None,
                      CNN_thres=0.0,
                      makecutout=True,
                      size=100,
                      size_cnn=32,
                      fmt='png',
                      outLevel=1,
                      nb_threads=8,
                      combined=False):
    """Filter transient candidates"""
    print('Filter candidates')
    # Take first candidate to extract the path where to store the file
    # No need to chack if substraction was performed,
    # as if it did only the ones from substracted files are with 'Match' == Y
    path, fname_ext = os.path.split(sources['filenames'][0])
    # Get rid of the extension to keep only the name
    fname2, extension = os.path.splitext(fname_ext)
    # Get rid of the _reg pattern
    fname2 = fname2.split('_ref')[0]

    # First get the sources not crossmatching with sources in catalogs
    mask_cat = sources['Match'] == 'N'

    # Remove candidates on the edges
    mask_edge = sources["edge"] == 'N'

    # Remove sources with FWHM ratio outside the desired range
    FWHM_ratio = sources["FWHM"] / sources["FWHMPSF"]
    mask_FWHM = (FWHM_ratio >= FWHM_ratio_lower) & \
                (FWHM_ratio <= FWHM_ratio_upper)

    mask_tot = mask_cat & mask_edge & mask_FWHM
    # Use a trained CNN model to filter candidates.
    if CNN_model is not None:
        print('Create fits cutouts for CNN')
        # Create fits cutouts to be be given to the CNN model
        path_CNN_cutouts = os.path.join(path, 'CNN_cutouts')
        mkdir_p(path_CNN_cutouts)
        args_data = []
        outnames = []
        info_dicts = []
        for cand in sources:
            coords = [cand['_RAJ2000'], cand['_DEJ2000']]
            outname = os.path.join(path_CNN_cutouts,
                                   'candidate_%d.fits' % (cand['idx']))
            info_dict = {}
            info_dict['RA'] = cand['_RAJ2000']
            info_dict['DEC'] = cand['_DEJ2000']
            info_dict['XPOS'] = cand['Xpos']
            info_dict['YPOS'] = cand['Ypos']
            info_dict['FILE'] = cand['filenames']
            info_dict['CANDID'] = cand['idx']
            info_dict['MAG'] = cand['mag_calib']
            info_dict['MAGERR'] = cand['mag_calib_err']
            info_dict['FWHM'] = cand['FWHM']
            info_dict['FWHMPSF'] = cand['FWHMPSF']

            args_data.append([cand['filenames'],
                              coords,
                              "world",
                              [size_cnn, size_cnn],
                              -1,
                              ])
            outnames.append(outname)
            info_dicts.append(info_dict)
        """
        # Run the make_sub_image in asynchroneous parallel
        # using many processes.
        # Need to cut the args list in the number of required threads
        # There is may be another way?
        N_sources = len(sources)
        # Check if there less data than number of threads.
        if N_sources < nb_threads:
            nb_threads = 1
        Ncut = int(N_sources / nb_threads)
        args_threads = []
        idx_stop = []
        for i in range(nb_threads):
            if i == 0:
                args_threads.append(args[i * Ncut : (i+1) * Ncut])
                idx_stop.append((i+1) * Ncut)
            elif i > 0 and i < nb_threads-1:
                args_threads.append(args[i * Ncut : (i+1) * Ncut])
                idx_stop.append((i+1) * Ncut)
            elif i == nb_threads-1:
                args_threads.append(args[i * Ncut : N_sources])
                idx_stop.append(N_sources)
            if idx_stop[-1] >= N_sources:
                break
        args_threads = np.array(args_threads)
        """
        args_data = np.array(args_data)
        pool = mp.Pool(nb_threads)
        # call apply_async() without callback
        """
        result_objects = [pool.apply_async(make_sub_image,
            args=(j[0,:], j[1,:], j[2,:], j[3,:], j[4,:]))
            for j in args_threads]
        """
        result_objects = [pool.apply_async(make_sub_image,
                                           args=(args_data[:, 0],
                                                 args_data[:, 1],
                                                 args_data[:, 2],
                                                 args_data[:, 3],
                                                 args_data[:, 4]))]

        # result_objects is a list of pool.ApplyResult objects
        results = [r.get() for r in result_objects]

        # Don't forget to close
        pool.close()
        pool.join()
        results = np.array(results[0])
        # Create fits cutouts
        p = mp.Pool(nb_threads)
        args = [[a, b, c, d, e, f] for a, b, c, d, e, f in zip(
            results[0, :],
            outnames,
            results[1, :],
            results[2, :],
            results[3, :],
            info_dicts)]
        p.starmap(make_fits, args)
        p.close()

        # Run CNN model to associate a probability to each cutout
        # The size of the cutout should be the same as the ones used
        # for the CNN training
        print("Use trained CNN model")
        infer(path_CNN_cutouts, CNN_model, 0.1)

        # Add the probability to the canidates table.
        infer_table = ascii.read(os.path.join(path_CNN_cutouts,
                                              'infer_results.dat'))
        sources = join(sources,
                       infer_table['idx', 'label0', 'label1'],
                       join_type='left')

        # keep only transients that are above the threshold
        mask_CNN = sources['label1'] >= CNN_thres
        mask_tot = mask_tot & mask_CNN

    # Write output file.
    candidates = sources[mask_tot]
    # Create ID to start from 1.
    candidates['cand_ID'] = np.arange(len(candidates)) + 1
    # Rename colums
    if 'label0' in candidates.colnames:
        candidates.rename_column('label0', 'P_False')
        candidates.rename_column('label1', 'P_True')
    candidates.write(
        os.path.join(path, fname2 + '_candidates.dat'),
        format='ascii.commented_header',
        overwrite=True
    )

    # Update
    print('Make cutouts')
    # Extract small image centered on candidates passing the filters
    if makecutout:
        path_cutout = os.path.join(path, 'cutouts')
        mkdir_p(path_cutout)
        args_data = []
        outnames = []
        if combined:
            args_combined = []
            path_cutout_combined = os.path.join(path_cutout,
                                                'combined')
            mkdir_p(path_cutout_combined)
        for cand in candidates:
            coords = [cand['_RAJ2000'], cand['_DEJ2000']]
            outname = os.path.join(
                path_cutout,
                'candidate_%d.%s' % (cand['cand_ID'], fmt)
            )
            if combined:
                outname_combined = os.path.join(
                    path_cutout_combined,
                    'candidate_%d_comb.%s' % (cand['cand_ID'], 'png'))

            header = fits.getheader(cand['OriginalIma'])
            try:
                date = Time(header["DATE-OBS"], format="fits")
                # convert in GPS time
                # date_JD = date.jd
            except BaseException:
                date = Time(header["MJD-OBS"], format="mjd")
            date.format = 'iso'
            if fmt != 'fits' or combined:
                _coords = SkyCoord(cand['_RAJ2000'], cand['_DEJ2000'],
                                   unit=(u.degree, u.degree),
                                   frame='icrs')
                coords_sexa = _coords.to_string(style='hmsdms')
                title = "RA Dec: %s \n" % coords_sexa + \
                        "Time (UTC): %s \n" % (date.value) + \
                        "Mag: %.2f +/- %.2f     " % (cand['mag_calib'],
                                                     cand['mag_calib_err']) + \
                        "     FWHM_ratio: %.2f" % (
                            cand['FWHM']/cand['FWHMPSF'])
                if CNN_model is not None:
                    title += "     CNN proba: %.2f " % cand['P_True']

            args_data.append([cand['filenames'],
                              coords,
                              "world",
                              [size_cnn, size_cnn],
                              -1,
                              ])
            outnames.append(outname)

            if combined:
                args_combined.append([
                    [cand['OriginalIma'],
                     cand['RefIma'],
                     cand['filenames']
                     ],
                    coords,
                    "world",
                    outname_combined,
                    [size, size],
                    -1,
                    title])

        # Create sub-array
        args_data = np.array(args_data)
        pool = mp.Pool(nb_threads)
        # call apply_async() without callback
        result_objects = [pool.apply_async(make_sub_image,
                                           args=(args_data[:, 0],
                                                 args_data[:, 1],
                                                 args_data[:, 2],
                                                 args_data[:, 3],
                                                 args_data[:, 4]))]

        # result_objects is a list of pool.ApplyResult objects
        results = [r.get() for r in result_objects]
        # Don't forget to close
        pool.close()
        pool.join()
        results = np.array(results[0])
        # Create fits cutouts
        p = mp.Pool(nb_threads)
        if fmt != 'fits':
            args = [[a, b, c, d, e] for a, b, c, d, e in zip(
                results[0, :],
                outnames,
                results[4, :],
                [fmt] * len(candidates),
                title)]
            p.starmap(make_figure, args)
        elif fmt == 'fits':
            args = [[a, b, c, d, e, f] for a, b, c, d, e, f in zip(
                results[0, :],
                outnames,
                results[1, :],
                results[2, :],
                results[3, :],
                info_dicts)]
            p.starmap(make_fits, args)
        p.close()

        if combined:
            print('Make combined cutouts')
            p = mp.Pool(nb_threads)
            p.starmap(combine_cutouts, args_combined)
            p.close()
