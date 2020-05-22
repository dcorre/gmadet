#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python-Pyraf module for GRANDMA detection of Optical candidates
Authors: David Corre, Orsay, France, corre@lal.in2p3.fr
          Martin Blazek, Granada, Spain, alf@iaa.es

Input arguments are filename, typical fwhm (or estimated by psfex), sextracting threshold and maximal distance for catalogue crosschecking in pixels
 
Example :
   python gmadet.py --filename /folder/image.fits --FWHM psfex --threshold 4 --radius_crossmatch 2.5 --telescope TRE 

"""


import sys, subprocess, glob, math, shutil, os
import argparse
import warnings

from catalogues import *
from phot_calibration import phot_calib
from utils import (load_config, clean_folder, send_data2DB, cut_image,
                   mv_p, mkdir_p, make_copy, clean_outputs)
from sanitise import sanitise_fits
from remove_cosmics import run_lacosmic
from astrometry import astrometrynet, scamp
from psfex import psfex
from substraction import substraction
from background import bkg_estimation

from astropy.io import ascii, fits
from astropy.table import vstack, Table, Column

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

from copy import deepcopy

warnings.simplefilter(action='ignore', category=FutureWarning)

def astrometric_calib(filenames, config, soft='scamp', accuracy=0.5, itermax=4, verbose='NORMAL'):
    """perform astrometric calibration"""

    imagelist = np.atleast_1d(filenames)
    for ima in imagelist:
        # Use scamp for astrometric calibration
        if soft == 'scamp':
            from astrometry import scamp
            scamp(ima, config, accuracy=accuracy, itermax=itermax, verbose=verbose)

        # Use astrometry.net for astrometric calibration
        elif soft == 'astrometrynet':
            from astrometry import astrometrynet

            # Get pixel scale in degrees
            header = fits.getheader(ima)
            try:
                pixScale = abs(header['CDELT1'])
            except Exception:
                try:
                    pixScale = abs(header['CD1_1'])
                except Exception:
                    print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1')
            # Set up boundaries for plate scale for astrometry.net
            scaleLow = 0.7 * pixScale * 3600
            scaleHigh = 1.3 * pixScale * 3600
            radius = max(header['NAXIS1']*pixScale, header['NAXIS2']*pixScale)
            asrometrynet(ima, radius=radius, scaleLow=scaleLow, scaleHigh=scaleHigh)


def sextractor(filelist, FWHM_list, thresh, telescope, config,
               verbose='NORMAL', subFiles=None, outLevel=1):
    """Run sextractor """
    # if substraction have been performed
    # Run sextractor on input image to get photometry calibrated and 
    # on substracted image to get interesting sources
    if subFiles:
        subFiles = np.array(subFiles)

        # No mask on iput data
        #mask = ['None'] * len(subFiles[:, 0])
        mask = ['None'] * len(filelist)
        mask.extend([im for im in subFiles[:, 3]])
        #weight_type = ['NONE'] * len(subFiles[:, 0])
        weight_type = ['NONE'] * len(filelist)
        weight_type.extend(['MAP_WEIGHT'] * len(mask))
       
        psfs = [im.split('.')[0] +'.psf' for im in filelist]
        # assume we take the psf of the original file
        psfs.extend([psfs[0]] * len(subFiles[:, 2]))
        #filelist = [im for im in subFiles[:, 0]]
        # Rather take the original before registration as it introduces artefact
        filelist = [im for im in filelist]
        
        filelist.extend([im for im in subFiles[:, 2]])
        # Duplicate FWHM list
        #FWHM_list.extend(FWHM_list)
        FWHM_list = np.ravel([[i]*len(filelist) for i in FWHM_list])
    else:
        mask = ['None'] * len(filelist)
        weight_type = ['NONE'] * len(filelist)
    for i, filename in enumerate(filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]
        if outLevel == 2:
            checkimage_type = 'BACKGROUND, SEGMENTATION'
            checkimage_name = folder + filename2 + '_background.fits' + ', ' + folder + filename2 + '_segmentation.fits'

        else:
            checkimage_type = 'NONE'
            checkimage_name = ' '

        subprocess.call(['sex', '-c', config['sextractor']['conf'], \
                filename, \
                '-WEIGHT_TYPE', str(weight_type[i]), \
                '-WEIGHT_IMAGE', str(mask[i]), \
                '-SEEING_FWHM', str(FWHM_list[i]), \
                '-DETECT_THRESH', str(thresh), \
                '-PARAMETERS_NAME', config['sextractor']['param'], \
                #'-FILTER_NAME', config['sextractor']['default_conv'], \
                '-CHECKIMAGE_TYPE', checkimage_type, \
                '-CHECKIMAGE_NAME',  checkimage_name, \
                '-VERBOSE_TYPE', verbose, \
                '-PSF_NAME', psfs[i], \
                '-CATALOG_NAME', folder + filename2 + '_SourcesDet.cat' ])


def get_photometry(filelist,FWHM_list,thresh,mkiraf=True,subFiles=None):
    """
    Performs sextracting by daofind and photometry by daophot 
    filename is WITHOUT suffix .fits
    FWHM_list is the estimation of FWHM for each image / quadrants
    THRESH is signal-to-noise ratio, typically 4 or so
    Outputs are two files *.coo.1 and *.mag.1 with x-y coordinates
    """
    
    # Create login.cl at execution of the script if flag set to true
    if mkiraf:
        proc = subprocess.Popen(['mkiraf'], stdin = subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outs, errs = proc.communicate('y\nxterm')

    # Now we make the Pyraf imports
    from pyraf import iraf
    from pyraf.iraf import daophot

    # Parameter to get rid of the faint stars with high magnitude error
    magnitude_error_threshold = 0.5

    iraf.noao
    iraf.digiphot
    iraf.daophot

    iraf.unlearn('phot')
    iraf.unlearn('datapars')
    iraf.unlearn('photpars')
    iraf.unlearn('daopars')

    iraf.datapars.epadu = 1
    #iraf.datapars.gain = 'GAIN'
    #iraf.datapars.ccdread = 'RON'
    iraf.datapars.readnoi = 6
    iraf.datapars.datamin = "INDEF"
    iraf.datapars.datamax = "50000"
    #iraf.datapars.datamax = "INDEF"


    for i, filename in enumerate(filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # fwhm tarot 1.5   oaj 3.5
        #iraf.datapars.fwhmpsf = FWHM_list[i]
        iraf.datapars.fwhm = FWHM_list[i]

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]
        daofind_output = folder + filename2 + '.coo.1'
        daophot_output = folder + filename2 + '.mag.1'

        IMstatres = iraf.imstat(filename,Stdout=1)
        IMmean = IMstatres[1].split()[2]
        iraf.datapars.sigma = math.sqrt(float(IMmean) * float(iraf.datapars.epadu) + float(iraf.datapars.readnoi)**2) / float(iraf.datapars.epadu)

        print("--- performing daophot sextracting ---")

        iraf.daofind(filename,output=daofind_output,verify="no", verbose="no", threshold=thresh)
        iraf.datapars.datamax = "INDEF"
        print("--- performing daophot photometry ---")
        iraf.daophot.phot(image=filename,coords=daofind_output, output=daophot_output, interactive="no", sigma="INDEF", airmass="AIRMASS", exposure="EXPOSURE", filter="FILTER", obstime="JD", calgorithm="gauss", verify="no", verbose="no")
    

def select_good_stars(filelist,limiting_mag_err,soft,edge_cut=32,sigma=1,subFiles=None):
# Performs selection of stars without INDEF in magnitude or magnitude error and with the flag "NoError" inside *.mag.1 file
# Saves into *.magfiltered file in x-y coordinates
# filename is WITHOUT suffix .fits
    # if substraction have been performed
    # Run sextractor on input image to get photometry calibrated and 
    # on substracted image to get interesting sources
    if subFiles:
        subFiles = np.array(subFiles)
        #filelist = [im for im in subFiles[:, 0]]
        # Rather take the original before registration as it introduces artefact
        original_filelist = list(filelist)
        original_filelist.extend([im for im in filelist]*len(subFiles))
        #originallist = [im for im in filelist]
        #filelist = []
        #for im, sub in zip(originallist, subFiles[:,2]):
        #    filelist.append(im)
        #    filelist.append(sub)
        subfiles = np.array(subFiles)
        filelist = [im for im in filelist]
        filelist.extend([im for im in subfiles[:, 2]])
    for filename in filelist:
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]
        fileext = '.' + filename_ext.split('.')[1]

        if soft == 'iraf':
            print ('Selecting only good stars from daophot output.\n')
            magfile = folder + filename2 + '.mag.1'
            resmaggile = folder + filename2 + ".magfiltered"
            #data = ascii.read(magfile)
            #data.show_in_browser()

            f1 = open(magfile, "r")
            f2 = open(resmaggile,"w")

            for kk in range(1,77):
                lajna = f1.readline()

            while lajna:
                lajna = f1.readline()
                xpos = lajna.split()[0]
                ypos = lajna.split()[1]
                lajna = f1.readline()
                lajna = f1.readline()
                lajna = f1.readline()
                mag = lajna.split()[4]
                merr = lajna.split()[5]
                if len(lajna.split()) < 8:
                    errmessage = lajna.split()[6]
                else:
                    errmessage = lajna.split()[7]
        
                if (errmessage != "NoError") or (mag == "INDEF") or (merr == "INDEF"):
                    #print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                    pass
                else:
                    #print(merr)
                    if (float(merr) < limiting_mag_err):
                        #print(xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                        f2.write(xpos+" "+ypos+" "+mag+" "+merr+"\n")                
                    else:
                        #print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                        pass

                lajna = f1.readline()

            f1.close()
            f2.close()

        elif soft == 'sextractor':
            sources = ascii.read(folder + filename2 + '_SourcesDet.cat', format='sextractor')
            mv_p(folder + filename2 + '_SourcesDet.cat',folder + filename2 + '_SourcesDetnoFilter.cat')
            # only if there is at least one detection
            if sources:
                # Remove sources too close to the imge edges
                header = fits.getheader(folder + filename2 + fileext)
                imsize = [int(header['NAXIS1']), int(header['NAXIS2'])]
                mask_edge = (sources['X_IMAGE'] > edge_cut) & \
                            (sources['Y_IMAGE'] > edge_cut) & \
                            (sources['X_IMAGE'] < imsize[1] - edge_cut) & \
                            (sources['Y_IMAGE'] < imsize[0] - edge_cut)

                """
                # Remove sources that are likely cosmic rays
                # Compute flux ratio as total_flux / nb_pixels for each source
                # Compute the median and std for sources with 10 < nb_pixels < 100
                # Higher than 10 to discard cosmics
                # Smaller than 100 to discard saturated sources
                # Remove sources with nb_pixels < 10 and fluxratio > fluxratio_med + fluxratio_std
                flux = sources['FLUX_AUTO']
                nbpix = sources['ISOAREA_IMAGE']
                fluxratio = flux / nbpix
                # Compute the median and std on the original image only
                # This assumes that the substracted image follows the original image in the for loop
                if '_sub' not in filename2:
                    mask = (nbpix > 10) & (nbpix < 100)
                    fluxratio_med = np.median(fluxratio[mask])
                    fluxratio_std = np.std(fluxratio[mask])

                mask_cosmics = (nbpix < 10) & (fluxratio > fluxratio_med + sigma * fluxratio_std)
                mask_tot = np.bitwise_and(mask_edge, np.invert(mask_cosmics))
                """

                # Remove sources too close to the edges
                #sources_filt = sources[mask_edge]
                # Flag sources too close to the edges
                # Only if there is at least one detection
                edge_flag = np.array(['N'] * len(sources))
                edge_flag[~mask_edge] = 'Y'
                sources['edge'] = edge_flag
            else:
                # Add something to make script not crashing
                sources['edge'] = []
            sources.write(folder + filename2 + '_SourcesDet.cat',format='ascii.commented_header', overwrite=True)


def convert_xy_radec(filelist, soft='sextractor', subFiles=None):
    """
    Performs pyraf transformation of x-y into RA-DEC coordinates
    filename is WITHOUT suffix .fits
    Input is the *.magfiltered file from select_good_stars() function
    Output is the *.magwcs file
    """
    """
        subFiles = np.array(subFiles)

        # No mask on iput data
        #mask = ['None'] * len(subFiles[:, 0])
        mask = ['None'] * len(filelist)
        print (len(subFiles[:, 0]))
        mask.extend([im for im in subFiles[:, 3]])
        #weight_type = ['NONE'] * len(subFiles[:, 0])
        weight_type = ['NONE'] * len(filelist)
        weight_type.extend(['MAP_WEIGHT'] * len(mask))

        #filelist = [im for im in subFiles[:, 0]]
        # Rather take the original before registration as it introduces artefact
        filelist = [im for im in filelist]
        filelist.extend([im for im in subFiles[:, 2]])
        # Duplicate FWHM list
        #FWHM_list.extend(FWHM_list)
        FWHM_list = np.ravel([[i]*len(filelist) for i in FWHM_list])
"""
    # If substraction has been performed
    if subFiles:
        original_filelist = list(filelist)
        original_filelist.extend([im for im in filelist]*len(subFiles))
        subfiles = np.array(subFiles)
        reference_filelist = ['None'] * len(filelist)
        reference_filelist.extend([im for im in subfiles[:, 1]])
        #filelist = [im for im in subfiles[:, 0]]
        # Rather take the original before registration as it introduces artefact
        filelist = [im for im in filelist]
        filelist.extend([im for im in subfiles[:, 2]])
    else:
        original_filelist = filelist
        reference_filelist = ['None']*len(filelist)

    for filename, original_filename, refimage  in zip(filelist, original_filelist, reference_filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]

        magfilewcs = folder + filename2 + ".magwcs"

        if soft == 'iraf':
            from pyraf import iraf

            magfile = folder + filename2 + ".magfiltered"
            iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")

            data1 = ascii.read(magfile, names=['Xpos', 'Ypos', 'Mag_aper', 'Mag_err_aper' ])
            #header = fits.getheader(filename)
            #w = wcs.WCS(header)
            #ra, dec = w.wcs_pix2world(data1['Xpos'], data1['Ypos'], 1)

            data2 = ascii.read(magfilewcs, names=['RA', 'DEC', 'Mag_aper', 'Mag_err_aper' ])
            data = Table([data1['Xpos'],data1['Ypos'], data2['RA'], data2['DEC'], data2['Mag_aper'], data2['Mag_err_aper'], [filename]*len(data1)], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_inst', 'Magerr_inst', 'filenames'])
            #data = Table([data1['Xpos'],data1['Ypos'], ra, dec, data1['Mag_aper'], data1['Mag_err_aper'], [filename]*len(data1)], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_inst', 'Magerr_inst', 'filenames'])
 
        elif soft == 'sextractor':
            sources = ascii.read(folder + filename2 + '_SourcesDet.cat')
            # If there is at least one detection
            if sources:
                header = fits.getheader(filename)
                w = wcs.WCS(header)
                ra, dec = w.wcs_pix2world(sources['X_IMAGE'], sources['Y_IMAGE'], 1)

                filenames = [filename] * len(ra)
            else:
                ra = []
                dec = []
                filenames = []
            data = Table([sources['X_IMAGE'],  sources['Y_IMAGE'], ra,dec, sources['MAG_AUTO'], sources['MAGERR_AUTO'], sources['edge'], sources['CHI2_PSF'], sources['CHI2_MODEL'], sources['FWHM_IMAGE'], sources['FWHMPSF_IMAGE'], sources['FLAGS_MODEL'], filenames], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_inst', 'Magerr_inst', 'edge', 'psf_chi2', 'model_chi2', 'FWHM', 'FWHMPSF', 'flag_psf', 'filenames'])
        
        # Flag to identify substraction image
        if '_sub' in magfilewcs:
            data['FlagSub'] = ['Y'] * len(data)
            data['OriginalIma'] = [original_filename]*len(data)
            data['RefIma'] = [refimage] * len(data)
        else:
            data['FlagSub'] = ['N'] * len(data)
            data['OriginalIma'] = [original_filename] * len(data)
            data['RefIma'] = [refimage] * len(data)
        
        
        data.write(magfilewcs, format='ascii.commented_header', overwrite=True)
        check_RADEC=data['RA', 'DEC']
        check_RADEC.write(magfilewcs+'2',format='ascii.commented_header', overwrite=True)


def crosscheck_with_catalogues(image_table, radius, catalogs=['I/345/gaia2', 'II/349/ps1', 'I/271/out', 'I/284/out'], Nb_cuts=(1,1), subFiles=None):
#def crosscheck_with_catalogues(image_table, radius, catalogs=['I/345/gaia2','II/349/ps1', 'I/271/out'], Nb_cuts=(1,1)):
    """
    Performs crosscheck with USNO B1.0 catalogue with *.magwcs
    filename is WITHOUT suffix .fits and maximal allowed difference radius is in arcseconds
    Input file is *.magwcs and the output is the list of the stars *.oc which were not identified in the catalogue
    radius is expressed in pixels

    Parameter for catalogue crosschecking, maximal allowed distance between sextracted and catalogue position in degrees
    Tarot 1px = 0.000907203 deg
    OAJ 1px = 1.543390792967E-04 deg
    AZT8 1px = 2.6321694198146E-04 deg
    allowed_crosscheck_radius = 3*0.000907203       # Tarot
    allowed_crosscheck_radius = 3*0.0001543390      # OAJ
    """

    cat_dict = {'I/284/out':'USNO-B1', 'I/345/gaia2':'GAIA DR2', 'II/349/ps1':'PS1 DR1', 'I/271/out':'GSC'}
    counter = 0
    _filename_list = []

    if subFiles is not None:
        subfiles = np.array(subFiles)
        #filelist = [im for im in subfiles[:, 0]]
        # Rather take the original data as the registartion introduces artefact
        filelist = [im for im in image_table['filenames']]
        filelist.extend([im for im in subfiles[:, 2]])
    else:
        filelist = image_table['filenames']
    for i, filename in enumerate(filelist):
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]
        extension = ''
        for ext in filename_ext.split('.')[1:]:
            extension = extension + '.' + ext

        magfilewcs = folder + filename2 + ".magwcs"

        if Nb_cuts == (1,1):
            original_name = folder + filename2
            quadrant = 1
        else:
            split_file = filename2.split('_Q')
            original_name = folder + split_file[0]
            for name in  split_file[1].split('_')[1:]:
                original_name = original_name + '_'  + name
            original_name = original_name.split('.')[0] 
            quadrant = split_file[1].split('_')[0]
            """
            if len(split_file[-1]) == 2:
                idx = 3
            elif len(split_file[-1]) == 3:
                idx = 4
            """
        _filename_list.append(original_name + ".oc")

        header = fits.getheader(filename)
        # Get pixel scale in degrees
        try:
            pixScale = abs(header['CDELT1'])
        except Exception:
            try:
                pixScale = abs(header['CD1_1'])
            except Exception:
                print ('Pixel scale could not be found in fits header.\n Expected keyword: CDELT1, _DELT1 or CD1_1')

        # Load detected sources in astropy table
        detected_sources = ascii.read(magfilewcs, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'edge', 'psf_chi2', 'model_chi2', 'FWHM', 'FWHMPSF', 'flag_psf', 'filenames', 'FlagSub', 'OriginalIma', 'RefIma'])
        if detected_sources:
            detected_sources['quadrant'] = [quadrant]*len(detected_sources)
            # Do not need it as the astrometric calibration is performed on each quadrant now.
            """
            # Transform X_pos and Y_pos to original image in case it was split
            header2 = fits.getheader(filename)
            Naxis1 = float(header2['NAXIS1'])
            Naxis2 = float(header2['NAXIS2'])
            Naxis11 = int(Naxis1/Nb_cuts[0])
            Naxis22 = int(Naxis2/Nb_cuts[1])
            if Nb_cuts == (1,1):
                quad = None
                index_i = 0
                index_j = 0
            else:
                quad, index_i, index_j = image_table['quadrant'][i].split('_')
                quad = quad[1:]
            detected_sources['Xpos'] = detected_sources['Xpos'] + Naxis22 * int(index_j)
            detected_sources['Ypos'] = detected_sources['Ypos'] + Naxis11 * int(index_i)
            """
            detected_sources['Xpos_quad'] = detected_sources['Xpos']
            detected_sources['Ypos_quad'] = detected_sources['Ypos']
        else:
            detected_sources['quadrant'] = []
            detected_sources['Xpos_quad'] = []
            detected_sources['Ypos_quad'] = []
        if i == 0:
            detected_sources_tot = deepcopy(detected_sources)
        else:
            detected_sources_tot = vstack([detected_sources_tot, detected_sources])
    # Add units
    detected_sources_tot['_RAJ2000'] *= u.deg
    detected_sources_tot['_DEJ2000'] *= u.deg
    # Add index for each source
    detected_sources_tot['idx'] = np.arange(len(detected_sources_tot))

    # Add a flag to indicate whether a source is crossmatched with 
    # a referenced source
    detected_sources_tot['Match'] = ['N'] * len(detected_sources_tot)

    # Initialise candidates with all detected sources
    candidates = deepcopy(detected_sources_tot)
    #candidates.write('test0.dat', format='ascii.commented_header', overwrite=True)
    print ('\nCrossmatching sources with catalogs.')
    print ('Radius used for crossmatching with catalogs: %.2f arcseconds\n' % (radius*pixScale*3600))

    mask_matched = candidates['Match'] == 'N'
    for catalog in catalogs:
        print (catalog, len(candidates[mask_matched]))
        # Use Xmatch to crossmatch with catalog
        crossmatch = run_xmatch(candidates[mask_matched], catalog, radius*pixScale*3600)

        #crossmatch.write('test.dat', format='ascii.commented_header', overwrite=True)
        # Do not consider duplicates
        # Meaning that if there are several sources
        # we consider it as a crossmatch
        referenced_star_idx = np.unique(crossmatch['idx'])

        candidates['Match'][referenced_star_idx ] = 'Y'

        # Update Match mask
        mask_matched = candidates['Match'] == 'N'

        if subFiles is not None:
            mask = (candidates['FlagSub'] == 'Y') & (mask_matched)
            mask2 = candidates['FlagSub'] == 'Y'
            print ('%d/%d candidates left in substracted image after crossmatching with %s' % (len(candidates[mask]),len(candidates[mask2]), cat_dict[catalog]))
            mask = (candidates['FlagSub'] == 'N') & (mask_matched)
            mask2 = np.invert(mask2)
            print ('%d/%d candidates left in original image (without substraction) after crossmatching with %s' % (len(candidates[mask]),len(candidates[mask2]), cat_dict[catalog]))
        else:
            print ('%d/%d candidates left after crossmatching with %s' % (len(candidates[mask_matched]),len(candidates), cat_dict[catalog]))
        if (len(candidates[mask_matched]) == 0):
            break

    # Get filename    
    _filename = np.unique(_filename_list)[0]
    # Write candidates file.
    # If substraction was performed, split transients into specific files
    if subfiles is not None:
        mask = (candidates['FlagSub'] == 'Y') & (mask_matched)
        candidates[mask].write(_filename.split('.')[0]+'_sub.oc', format='ascii.commented_header', overwrite=True)
    
    oc = candidates[mask]['_RAJ2000', '_DEJ2000']
    oc.write(_filename.split('.')[0]+'_sub.oc_RADEC',format='ascii.commented_header', overwrite=True)
    mask = (candidates['FlagSub'] == 'N') & (mask_matched)
    candidates[mask].write(_filename, format='ascii.commented_header', overwrite=True)
    oc=candidates[mask]['_RAJ2000', '_DEJ2000']
    oc.write(_filename.split('.')[0]+'.oc_RADEC',format='ascii.commented_header', overwrite=True)

    # Also write a file with all the sources detected to know
    candidates.write(_filename.split('.')[0] + '.alldetections', format='ascii.commented_header', overwrite=True)


    #oc=candidates['Xpos', 'Ypos']
    #oc.write(magfilewcs+'4',format='ascii.commented_header', overwrite=True)

    return candidates

def check_moving_objects(filelist):
    """
    Crossmatch the list of candidates with moving objects using SkyBoT
    
    """
    for filename in filelist:
        header = fits.getheader(filename)
        ra_deg = float(header['CRVAL1'])*u.deg
        dec_deg = float(header['CRVAL2'])*u.deg
        date = header['DATE-GPS']
        radius = 2 *u.deg
        Texp = float(header['exposure']) * u.second
    
        moving_objects = skybot(ra_deg, dec_deg, date, radius, Texp)

        moving_objects.write('moving_objects.dat', format='ascii.commented_header', overwrite=True)
        moving_objects['RA', 'DEC'].write('moving_objects.reg', format='ascii.commented_header', overwrite=True)

        candidates_out = crossmatch_skybot(candidates, moving_objects, radius=10)

        print ('%d match with a moving object found' % (len(candidates)-len(candidates_out)))

    #return candidates_out

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Finding unknown objects in astronomical images.')

    parser.add_argument('--filename',
                        dest='filename',
                        required=True,
                        type=str,
                        help='Path to file')

    parser.add_argument('--FWHM',
                        dest='FWHM',
                        required=True,
                        help='Typical telescope FWHM')

    parser.add_argument('--mag_err_cut',
                        dest='mag_err_cut',
                        required=False,
                        default=0.5,
                        type=float,
                        help='Consider only sources with magnitude error < mag_err_cut')

    parser.add_argument('--radius_crossmatch',
                        dest='radius_crossmatch',
                        required=True,
                        type=float,
                        help='Radius to use for crossmatching, in arcseconds')

    parser.add_argument('--threshold',
                        dest='threshold',
                        required=True,
                        type=float,
                        help='Consider only sources above this threshold')


    parser.add_argument('--soft',
                        dest='soft',
                        required=False,
                        choices=['sextractor', 'iraf'],
                        default='sextractor',
                        type=str,
                        help='Soft to use for detecting sources')


    parser.add_argument('--telescope',
                        dest='telescope',
                        choices=['TRE','TCA','TCH','OAJ-T80','Lisniky-AZT8','UBAI-T60S','UBAI-T60N', 'FRAM-CTA-N', 'FRAM-Auger', 'KAIT', 'IRIS'],
                        required=True,
                        type=str,
                        help='Alias for the telescopes')

    parser.add_argument('--owncloud_path',
                        dest='owncloud_path',
                        required=False,
                        type=str,
                        help='Local path to the owncloud')

    parser.add_argument('--VOE_path',
                        dest='VOE_path',
                        required=False,
                        type=str,
                        help='Path + filename of the VoEvent containing the observation plan.')

    parser.add_argument('--quadrants',
                        dest='quadrants',
                        required=False,
                        default=1,
                        type=int,
                        help='Number of quadrants the image is divided.')

    parser.add_argument('--doAstrometry',
                        dest='doAstrometry',
                        required=False,
                        default='scamp',
                        choices=['No', 'scamp', 'astrometrynet'],
                        type=str,
                        help='Whether to perform astrometric calibration, with scamp or astrometry.net.')

    parser.add_argument('--verbose',
                        dest='verbose',
                        required=False,
                        default='NORMAL',
                        choices=['QUIET', 'NORMAL', 'FULL', 'LOG'],
                        type=str,
                        help='Level of verbose, according to astromatic software')

    parser.add_argument('--doSub',
                        dest='doSub',
                        required=False,
                        type=str,
                        help='Whether to perform astrometric calibration, with ps1 images or user provided reference image. Type "ps1" for PS1 reference image or provide the path to your reference image.')

    parser.add_argument('--ps1_method',
                        dest='ps1_method',
                        required=False,
                        default='mosaic',
                        choices=['mosaic', 'individual'],
                        type=str,
                        help='When substracting images using Pan-STARRS reference images, there 2 options, either create a mosaic of all PS1 image and substract or do the substraction individually for each PS1 image. In the latter case, your image is cut to match the PS1 image.')


    parser.add_argument('--Remove_cosmics',
                        dest='Remove_cosmics',
                        action='store_true',
                        help='Whether to remove cosmic rays using lacosmic.')

    parser.add_argument('--sub_bkg',
                        dest='sub_bkg',
                        action='store_true',
                        help='Whether to substract background.')

    parser.add_argument('--output_data_level',
                        dest='outLevel',
                        required=False,
                        type=int,
                        default=1,
                        choices=[0,1,2],
                        help='Number of output files that are kept after the process.')

    args = parser.parse_args()

    Nb_cuts = (args.quadrants,args.quadrants)
    
    # Load config files for a given telescope
    config = load_config(args.telescope)

    # List all the files in the given path
    if os.path.isdir(args.filename):
        # expected extensions: .fits and .fit
        # Get all the prefixes corresponding to one field
        filenames = glob.glob(args.filename + '/**/*.fit*', recursive=True)
    else:
        filenames = [args.filename]

    # copy original images
    # Create list of the copy images
    filenames = make_copy(filenames, outputDir='gmadet_results/')



    for filename in filenames:

        print ('Sanitise header and data of %s.\n' % filename)
        sanitise_fits(filename)

        # Cut image into several quadrants if required
        # And create table with filename and quadrant ID
        image_table = cut_image(filename, config, Nb_cuts=Nb_cuts, doAstrometry=args.doAstrometry)

        if args.FWHM == 'psfex':
            # Estimate the PSF FWHM for each image/quadrants using psfex
            FWHM_list = psfex(image_table['filenames'], config, verbose=args.verbose, outLevel=args.outLevel)
        else:
            FWHM_list = [args.FWHM] * len(image_table)

        if args.sub_bkg:
            # Substract background
            bkg_estimation(image_table['filenames'], box=(20,20), filter_size=(3,3), outLevel=args.outLevel)

        if args.Remove_cosmics:
            print ('Running lacosmic on %s to remove cosmic rays. \n' % filename )
            # Clean cosmic rays
            run_lacosmic(image_table['filenames'], FWHM_list, sigma=5, outLevel=args.outLevel)

        if args.doAstrometry != 'No':
            astrometric_calib(image_table['filenames'], config, soft=args.doAstrometry, verbose=args.verbose, accuracy=0.5)

        if args.doSub:
            substracted_files = substraction(image_table['filenames'], args.doSub, config, soft='hotpants', method=args.ps1_method, verbose=args.verbose, outLevel=args.outLevel)
        else:
            substracted_files = None

        if args.soft == 'iraf':
            clean_folder(image_table['filenames'], subFiles=substracted_files)
            get_photometry(image_table['filenames'], FWHM_list, args.threshold, subFiles=substracted_files)
        elif args.soft == 'sextractor':
            sextractor(image_table['filenames'], FWHM_list, args.threshold, args.telescope, config,verbose=args.verbose, subFiles=substracted_files, outLevel=args.outLevel)

        select_good_stars(image_table['filenames'], args.mag_err_cut, args.soft, sigma=1, subFiles=substracted_files)
        convert_xy_radec(image_table['filenames'],soft=args.soft, subFiles=substracted_files)
        total_candidates = crosscheck_with_catalogues(image_table,args.radius_crossmatch, Nb_cuts=Nb_cuts, subFiles=substracted_files)
        #check_moving_objects(args.filename, total_candidates)

        #total_candidates = ascii.read('total_candidates.dat')
        total_candidates_calib = phot_calib(total_candidates, args.telescope, radius=args.radius_crossmatch,doPlot=True, subFiles=substracted_files)

        #total_candidates_calib = ascii.read('Test_sendDB/gmadet_results/jul1919-010r_sh_tot_cand2.dat')

        # If both arguments VOE_path and owncloud_path are provided
        # Send candidates to database
        # Set the tile_id corresponding to your tile by hand at the moment
        if args.VOE_path and args.owncloud_path:
            send_data2DB(filename, total_candidates_calib, Nb_cuts, args.owncloud_path, args.VOE_path, "utilsDB/usrpwd.json",debug=True, subFiles=substracted_files)

        # clean output files
        clean_outputs(image_table['filenames'], args.outLevel)

