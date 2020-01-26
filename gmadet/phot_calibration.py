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
from utils import get_filter

from astropy.io import ascii, fits
from astropy import units as u
from astropy.table import vstack, Table, Column

from copy import deepcopy

warnings.simplefilter(action='ignore')

def phot_calib(candidates_list, telescope, catalog='II/349/ps1', radius = 3, doPlot=True):
    """Perform photometric calibration using catalogs"""

    delta_mag_median_list = []
    filename_list = []

    # Get sources 
    for i, key in enumerate(candidates_list.group_by('filenames').groups.keys) :
        print ('Processing photometric calibration for ', key[0])

        # Get pixel scale in degrees
        header = fits.getheader(key[0])
        try:
            pixScale = abs(header['CDELT1'])
        except Exception:
            try:
                pixScale = abs(header['_DELT1'])
            except Exception:
                try:
                    pixScale = abs(header['CD1_1'])
                except Exception:
                    print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1')


        # Get filter
        #band_DB, band_cat = get_filter(key[0], telescope)
        band_DB = 'Clear'
        band_cat = 'g+r'

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
        fname = folder + fname2 + '.magwcs'

        detected_sources = ascii.read(fname, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames'])
        # Add units
        detected_sources['_RAJ2000'] *= u.deg
        detected_sources['_DEJ2000'] *= u.deg
        # Add index for each source
        detected_sources['idx'] = np.arange(len(detected_sources))

        crossmatch = run_xmatch(detected_sources, catalog, radius*pixScale*3600)
        # Initialise flag array. 0: unknown sources / 1: known sources
        flag = np.zeros(len(crossmatch))
        # Do not consider duplicates
        referenced_star_idx = np.unique(crossmatch['idx'])

        crossmatch['id'] = np.arange(len(crossmatch))
        closest_id = []
        for idx in referenced_star_idx:
            mask = crossmatch['idx'] == idx
            closest_id.append(crossmatch[mask]['id'][0])
        
        # Set flag indexes to 1 for detected sources associated to a star
        flag[closest_id] = 1


        ref_sources = crossmatch[flag == 1]

        bands = band_cat.split('+')
        # Bunch of conditions to filter the catalog data
        mask = ref_sources['Nd'] > 15
        ref_sources=ref_sources[mask]
        for band in bands:
            mask = ref_sources['%sFlags' % band] == 115000
            ref_sources=ref_sources[mask]

        #mask = ((ref_sources['Nd'] > 15) & (ref_sources['rFlags'] == 115000) & (ref_sources['gFlags'] == 115000))

        # Combine filter bands if needed 
        if len (bands) > 1:
            #jansky = []
            jansky = np.zeros(len(ref_sources))
            for band in bands:
                jansky = jansky + 3631 * 10**(-0.4*(ref_sources['%smag' % band]))
            newmag = -2.5*np.log10(jansky/(3631))
        else:
            newmag = ref_sources['%smag' % band]


        # Remove objects to get all objects with delte_mag < 1 sigma
        delta_mag = ref_sources['mag_inst'] - newmag
        delta_mag_median = np.median(delta_mag)
        delta_mag_std = np.std(delta_mag)
        counter = 1
        while (np.max(abs(delta_mag)) > abs(delta_mag_median) + delta_mag_std) and (counter <= 5):
            
            #print (len(delta_mag), delta_mag_median, delta_mag_std)
            mask = abs(delta_mag) < abs(delta_mag_median) + delta_mag_std
            ref_sources = ref_sources[mask]
            newmag = newmag[mask]
             
            delta_mag = ref_sources['mag_inst'] - newmag
            delta_mag_median = np.median(delta_mag)
            delta_mag_std = np.std(delta_mag)
            counter += 1

        #delta_mag_std = np.std(delta_mag)
        ref_sources.write(folder+fname2+'_ZP_%d.dat' % i, format='ascii.commented_header', overwrite=True)
        if doPlot:
            ref_sources.show_in_browser(jsviewer=True)
            plt.scatter(ref_sources['mag_inst'], ref_sources['rmag'], color ='blue')
            plt.scatter(ref_sources['mag_inst'], newmag, color='red')
            plt.plot(ref_sources['mag_inst'], ref_sources['mag_inst']-delta_mag_median, color='green')
            plt.savefig(folder+fname2+'_ZP_%d.png' % i)
            plt.show()

        delta_mag_median_list.append(delta_mag_median)
        filename_list.append(key[0])

    
    # Apply photmetric calibration to candidates
    mag_calib_col = Column(np.zeros(len(candidates_list)), name='mag_calib')
    mag_calib_err_col = Column(np.zeros(len(candidates_list)), name='mag_calib_err')
    magsys_col = Column(['None']*len(candidates_list), name='magsys')
    filter_cat_col = Column(['None']*len(candidates_list), name='filter_cat')
    filter_DB_col = Column(['NoFilterFound']*len(candidates_list), name='filter_DB')

    candidates_list.add_columns([mag_calib_col, mag_calib_err_col, magsys_col, filter_cat_col, filter_DB_col])

    for j, filename in enumerate(filename_list):
        mask = candidates_list['filenames'] == filename
        candidates_list['mag_calib'][mask] = candidates_list['mag_inst'][mask] - delta_mag_median
        # Quadratic sum of statistics and calibration errors. 
        candidates_list['mag_calib_err'][mask] = np.sqrt(candidates_list['mag_inst_err'][mask]**2 + delta_mag_std**2)

        # Define magnitude system
        if catalog == 'II/349/ps1':
            candidates_list['magsys'][mask] = 'AB'
        else:
            pass

        candidates_list['filter_cat'][mask] = band_cat
        candidates_list['filter_DB'][mask] = band_DB


    candidates_list.write(folder+fname2+'_tot_cand2.dat', format='ascii.commented_header', overwrite=True)


    return  candidates_list
