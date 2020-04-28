#! usr/bin/env python
# -*- coding: utf-8 -*-

# Python module for performing photometric calibration
# David Corre, corre@lal.in2p3.fr
#

import os
import numpy as np
import matplotlib.pyplot as plt
import warnings 

from catalogues import run_xmatch
from utils import get_phot_cat, filter_catalog_data
from phot_conversion import *

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import vstack, Table, Column
from astropy.stats import sigma_clip

from copy import deepcopy

warnings.simplefilter(action='ignore')

def crossmatch(detected_sources, radius, pixScale, catalog):
    """
    Crossmatch the output file of Sextractor / pyRAF
    Remove duplicate source during crossmatch
    Keep only the closest match
    """

    # Import Sextractor or pyRAF results
    #detected_sources = ascii.read(fname, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'edge', 'filenames', 'FlagSub', 'OriginalIma', 'RefIma'])

    # Add units
    detected_sources['_RAJ2000'] *= u.deg
    detected_sources['_DEJ2000'] *= u.deg
    # Add index for each source
    detected_sources['idx'] = np.arange(len(detected_sources))

    # Do not consider sources close to the image edges
    mask = detected_sources['edge'] == 'N'

    # Run xmatch to crossmatch detected sources with available catalog
    crossmatch = run_xmatch(detected_sources[mask], catalog, radius*pixScale*3600)
    # Initialise flag array. 0: unknown sources / 1: known sources
    flag = np.zeros(len(crossmatch))
    # Do not consider duplicates
    referenced_star_idx = np.unique(crossmatch['idx'])
    # Consider only closest crossmatch if multiple association.
    crossmatch['id'] = np.arange(len(crossmatch))
    closest_id = []
    for idx in referenced_star_idx:
        mask = crossmatch['idx'] == idx
        closest_id.append(crossmatch[mask]['id'][0])
    # Set flag indexes to 1 for detected sources associated to a star
    flag[closest_id] = 1
    # ref_sources is an astropy table containing a single crossmatch
    # for detected sources. Sources not crossmatched are not included.
    # It is used for photometric calibration purpose.
    ref_sources = crossmatch[flag == 1]

    return ref_sources

def conv_mag_sys(data, band, catalog):
    """
    For a gieven telescope magnitude and refence star catalog, perform
    conversion from catalog bands to telescope bands
    """
    if band in ['g', 'r', 'i', 'z', 'y', 'g+r']:
        if catalog == 'II/349/ps1':
            # No transformation
            bands = band.split('+')
            if len (bands) > 1:
                #jansky = []
                jansky = np.zeros(len(data))
                mag_err = np.zeros(len(data))
                for filt in bands:
                    jansky = jansky + 3631 * 10**(-0.4*(data['%smag' % filt]))
                    # Add error in quadrature. Might be too simple
                    mag_err = mag_err + data['e_%smag' % filt]**2
                newmag = -2.5*np.log10(jansky/(3631))
                newmag_err = np.sqrt(mag_err)
            else:
                newmag = data['%smag' % bands[0]]
                newmag_err = data['e_%smag' % bands[0]]

            data['mag_cat'] = newmag
            data['magerr_cat'] = newmag_err
            data['magsys'] = 'AB'
            catalogName = 'Pan-STARRS DR1'

        elif catalog == 'V/147/sdss12':
            # No transformation
            catalogName = 'SDSS DR12'
            pass
        elif catalog == 'I/345/gaia2':
            bands = band.split('+')
            if len (bands) > 1:
                #jansky = []
                jansky = np.zeros(len(data))
                mag_err = np.zeros(len(data))
                for filt in bands:
                    data2 = data = gaia2SDSS(filt, data)
                    jansky = jansky + 3631 * 10**(-0.4*(data2['%s_SDSSMag' % filt]))
                    # Add error in quadrature. Might be too simple
                    mag_err = mag_err + data2['calib_err']**2
                newmag = -2.5*np.log10(jansky/(3631))
                newmag_err = np.sqrt(mag_err)
                data['mag_cat'] = newmag
                data['magerr_cat'] = newmag_err
                data['magsys'] = 'AB'
            else:
                data = gaia2SDSS(band, data)
                data.rename_column('%s_SDSSMag' % band, 'mag_cat')
                data.rename_column('calib_err', 'magerr_cat')
                data['magsys'] = 'AB'
            catalogName = 'Gaia DR2'

    elif band in ['B', 'V', 'R', 'I']:
        if catalog == 'II/349/ps1':
            data = PS2Johnson(band, data)
            data.rename_column('%sMag' % band, 'mag_cat')
            data.rename_column('calib_err', 'magerr_cat')
            data['magsys'] = 'AB'
            catalogName = 'Pan-STARRS DR1'
        elif catalog == 'V/147/sdss12':
            data = SDSS2Johnson(band, data)
            data.rename_column('%sMag' % band, 'mag_cat')
            data.rename_column('calib_err', 'magerr_cat')
            data['magsys'] = 'AB'
            catalogName = 'SDSS DR12'
        elif catalog == 'I/345/gaia2':
            data = gaia2Johnson(band, data)
            data.rename_column('%sMag' % band, 'mag_cat')
            data.rename_column('calib_err', 'magerr_cat')
            data['magsys'] = 'Vega'
            catalogName = 'Gaia DR2'
        elif catalog == 'I/284/out':
            # Need to add
            catalogName = 'USNO-B1'
            pass

    return data, catalogName

def zeropoint(data, sigma, quadrant, folder, fname, band, catalog, doPlot=False):
    """"Compute zeropoints"""
    #data.show_in_browser()
    # Sigma clipping for zeropoints
    #clip = sigma_clip(a['Delta_Mag'], sigma = 1.5, masked=True)
    #clip_mask = np.invert(clip.recordmask)
    #a = a[clip_mask]
    # Remove objects to get all objects with delte_mag < 1 sigma
    delta_mag = data['mag_inst'] - data['mag_cat']

    clip = sigma_clip(delta_mag, sigma = sigma)
    clip_mask = np.invert(clip.recordmask)

    newdata = data[clip_mask]
    delta_mag = newdata['mag_inst'] - newdata['mag_cat']
    delta_mag_median = np.median(delta_mag)
    delta_mag_std = np.std(delta_mag)
    
    newdata.write(folder+fname+'_ZP_%d.dat' % quadrant, format='ascii.commented_header', overwrite=True)
    
    if doPlot:
        #newdata.show_in_browser(jsviewer=True)
        plt.figure()
        plt.scatter(data['mag_inst'], data['mag_cat'], color ='blue')
        median_nosigmaclip = np.median(data['mag_inst'] - data['mag_cat'])
        std_nosigmaclip = np.std(data['mag_inst'] - data['mag_cat'])
        plt.plot(data['mag_inst'], data['mag_inst']-median_nosigmaclip, color='green')
        plt.xlabel('Instrumental magnitude')
        plt.ylabel('%s band %s' % (band,catalog))
        plt.title('zeropoints without sigma clipping.\nMedian: %.2f, std: %.2f' % (median_nosigmaclip, std_nosigmaclip))
        plt.savefig(folder+fname+'_ZP_%d_nosigmaClipping.png' % quadrant)

        plt.figure()
        plt.scatter(newdata['mag_inst'], newdata['mag_cat'], color ='blue')
        plt.plot(newdata['mag_inst'], newdata['mag_inst']-delta_mag_median, color='green')
        plt.xlabel('Instrumental magnitude')
        plt.ylabel('%s band %s' % (band,catalog))
        plt.title('zeropoints after sigma clipping (sigma=1.5).\nMedian: %.2f, std: %.2f' % (-delta_mag_median,delta_mag_std))
        plt.savefig(folder+fname+'_ZP_%d.png' % quadrant)
        #plt.show()

    return newdata, delta_mag_median, delta_mag_std 

def phot_calib(candidates_list, telescope, radius = 3, doPlot=True, subFiles=None):
    """Perform photometric calibration using catalogs"""

    deltaMagMedianlist = []
    deltaMagStdlist = []
    filename_list = []

    # Filter to perfrom the calibration on the original image.
    # Not the substracted image
    #mask = candidates_list['FlagSub'] == 'N'

    #candidates_list.group_by('filenames').show_in_browser()
    #candidates_list[mask].group_by('filenames').show_in_browser()
    # Get sources 
    for i, key in enumerate(candidates_list.group_by('OriginalIma').groups.keys) :
        print ('Processing photometric calibration for ', key[0])

        # Get path and filename to images
        path, fname_ext = os.path.split(key[0])
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        fname2 = fname_ext.split('.')[0]
        extension = ''
        for ext in fname_ext.split('.')[1:]:
            extension = extension + '.' + ext

        # Get pixel scale in degrees
        header = fits.getheader(key[0])
        try:
            pixScale = abs(header['CDELT1'])
        except Exception:
            try:
                pixScale = abs(header['CD1_1'])
            except Exception:
                print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1 or CD1_1')

        # Get filter and catalog to perform photometric calibration
        band_DB, band_cat, catalog = get_phot_cat(key[0], telescope)

        # Import Sextractor or pyRAF results
        #fname = folder + fname2 + '.magwcs'
        fname = folder + fname2 + '.alldetections'
        detected_sources = ascii.read(fname)
        print ('Crossmatching with catalog.')
        ref_sources = crossmatch(detected_sources, radius, pixScale, catalog)

        # Remove extended sources and bad measurements from reference
        # stars catalog
        # ref_sources.show_in_browser(jsviewer=True)
        print ('Removed extended objects or with band quality flags.')
        good_ref_sources = filter_catalog_data(ref_sources, catalog)
        #good_ref_sources.show_in_browser(jsviewer=True)

        # Transform filter bands in catalog to telescope ones
        print ('Convert magnitude system.')
        ref_cat, catalogName = conv_mag_sys(good_ref_sources, band_cat, catalog)

        print ('Compute zeropoint.')
        ref_cat_calibrated, deltaMagMedian, deltaMagStd = zeropoint(ref_cat, 1.5, i, folder, fname2, band_cat, catalogName, doPlot=True)

        deltaMagMedianlist.append(deltaMagMedian)
        deltaMagStdlist.append(deltaMagStd)
        filename_list.append(key[0])

        # Adding calibrated mag to all detected sources
        detected_sources['ZP'] = [deltaMagMedian] * len(detected_sources)
        detected_sources['ZP_err'] = [deltaMagStd] * len(detected_sources)
        detected_sources['mag_calib'] = detected_sources['mag_inst'] - deltaMagMedian
        detected_sources['mag_calib_err'] = np.sqrt(detected_sources['mag_inst_err']**2 + deltaMagStd**2)
        detected_sources['filter_cat'] = [band_cat] * len(detected_sources)
        detected_sources['filter_DB'] = [band_DB] * len(detected_sources)

        detected_sources.write(fname, format='ascii.commented_header', overwrite=True)

    print ('Add magnitude to candidates.') 
    # Apply photmetric calibration to candidates
    mag_calib_col = Column(np.zeros(len(candidates_list)), name='mag_calib')
    mag_calib_err_col = Column(np.zeros(len(candidates_list)), name='mag_calib_err')
    magsys_col = Column(['None']*len(candidates_list), name='magsys')
    filter_cat_col = Column(['None']*len(candidates_list), name='filter_cat')
    filter_DB_col = Column(['NoFilterFound']*len(candidates_list), name='filter_DB')

    candidates_list.add_columns([mag_calib_col, mag_calib_err_col, magsys_col, filter_cat_col, filter_DB_col])

    for j, filename in enumerate(filename_list):
        # Compute calibrated magnitudes for transient candidates only
        mask = candidates_list['OriginalIma'] == filename
        candidates_list['mag_calib'][mask] = candidates_list['mag_inst'][mask] - deltaMagMedianlist[j]
        # Quadratic sum of statistics and calibration errors. 
        candidates_list['mag_calib_err'][mask] = np.sqrt(candidates_list['mag_inst_err'][mask]**2 + deltaMagStdlist[j]**2)

        # Define magnitude system
        if catalog == 'II/349/ps1':
            candidates_list['magsys'][mask] = 'AB'
        else:
            pass

        candidates_list['filter_cat'][mask] = band_cat
        candidates_list['filter_DB'][mask] = band_DB

    candidates_list.write(folder+fname2+'_transient_candidates.dat', format='ascii.commented_header', overwrite=True)

    return  candidates_list
