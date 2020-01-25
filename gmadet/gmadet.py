#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Python-Pyraf module for GRANDMA detection of Optical candidates
# Authors: David Corre, Orsay, France, corre@lal.in2p3.fr
#          Martin Blazek, Granada, Spain, alf@iaa.es
#
# Input arguments are filename, typical fwhm, sextracting threshold and maximal distance for catalogue crosschecking in degrees
# 
# Example:
#   python gmadet.py --filename /folder/image.fits --FWHM 3.5 --threshold 4 --radius_crossmatch 3 --soft iraf --telescope TRE 

import sys, subprocess, glob, math, shutil, os
import argparse
import warnings

from catalogues import *
from phot_calibration import phot_calib
from utils import load_config, clean_folder, send_data2DB, mv_p, mkdir_p
from astrometry import astrometrynet, scamp

from astropy.io import ascii, fits
from astropy.table import vstack, Table, Column

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

from copy import deepcopy

warnings.simplefilter(action='ignore', category=FutureWarning)


def astrometric_calib(filename, config, Nb_cuts=(2,2), soft='scamp', outputDir='gmadet_results/', accuracy=0.6, itermax=3, verbose='NORMAL'):
    """perform astrometric calibration"""
    from astrometry import astrometrynet, scamp
    from utils import cut_image 

    if soft in ['scamp', 'astrometrynet']:
        doAstrometry = True
    else:
        doAstrometry = False

    path, filename_ext = os.path.split(filename)

    # Get rid of the extension to keep only the name
    filename2 = filename_ext.split('.')[0]
    extension = ''
    for ext in filename_ext.split('.')[1:]:
        extension = extension + '.' + ext

    if path:
        folder = path + '/'
    else:
        folder = ''

    resultDir = folder+outputDir
    # Create results folder
    mkdir_p(resultDir)

    # Cut the image in the required number of quadrants
    cut_image(filename, resultDir, Nb_cuts = Nb_cuts)

    filelist = []
    quadrant = []
    
    if doAstrometry:
        # Use scamp for astrometric calibration
        if soft == 'scamp':

            if Nb_cuts[0] > 1 or Nb_cuts[1] > 1:
                index = 0
                for i in range(Nb_cuts[0]):
                    for j in range(Nb_cuts[1]):
                        index+=1
                        filein = resultDir + filename2 + "_Q%d" % (index) + extension
                        filelist.append(filein)
                        ii=0
                        doffset = 10
                        while (doffset >= accuracy) and (ii <= itermax):
                            ii+=1
                            doffset = scamp(filein, config, verbose=verbose)
                            print ('Astrometric precision run %d: %.2f arcseconds' % (ii, doffset))

                        quadrant.append('Q%d_%d_%d' % (index,i,j))
            else:
                filelist.append(resultDir + filename2 + extension)
                ii=0
                doffset = 10
                while (doffset >= accuracy) and (ii <= itermax):
                    ii+=1
                    doffset = scamp(resultDir + filename2 + extension, config, verbose=verbose)
                    print ('Astrometric precision run %d: %.2f arcseconds' % (ii, doffset))

                quadrant.append('None')

        # Use astrometry.net for astrometric calibration
        elif soft == 'astrometrynet':

            # Get pixel scale in degrees
            header = fits.getheader(filename)
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
            # Set up boundaries for plate scale for astrometry.net
            scaleLow = 0.7 * pixScale * 3600
            scaleHigh = 1.3 * pixScale * 3600
            radius = max(header['NAXIS1']*pixScale, header['NAXIS2']*pixScale)
            if Nb_cuts[0] > 1 or Nb_cuts[1] > 1:
                index = 0
                for i in range(Nb_cuts[0]):
                    for j in range(Nb_cuts[1]):
                        index+=1
                        filein = resultDir + filename2 + "_Q%d" % (index) + extension
                        filelist.append(filein)
                        asrometrynet(filein,radius=radius,scaleLow=scaleLow,scaleHigh=scaleHigh)
                        quadrant.append('Q%d_%d_%d' % (index,i,j))
            else:
                filelist.append(resultDir + filename2 + extension)
                astrometrynet(resultDir + filename2 + extension,radius=radius,scaleLow=scaleLow,scaleHigh=scaleHigh)
                quadrant.append('None')
    else:
        filelist.append(resultDir + filename2 + extension)
        quadrant.append('None')

    image_table = Table([filelist, quadrant], names=['filenames', 'quadrant']) 
    return image_table


def sextractor(filelist, fwhmpsf, thresh, telescope, config, verbose='NORMAL'):
    """Run sextractor """

    for filename in filelist:
        print ('sex',filename)
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]

        subprocess.call(['sex', '-c', config['sextractor']['conf'], \
                filename, \
                '-SEEING_FWHM', str(fwhmpsf), \
                '-DETECT_THRESH', str(thresh), \
                '-PARAMETERS_NAME', config['sextractor']['param'], \
                '-FILTER_NAME', config['sextractor']['default_conv'], \
                '-CHECKIMAGE_TYPE', 'SEGMENTATION', \
                '-CHECKIMAGE_NAME', folder + filename2 + '_segmentation.fits', \
                '-VERBOSE_TYPE', verbose, \
                '-CATALOG_NAME', folder + filename2 + '_SourcesDet.cat' ])



def get_photometry(filelist,fwhmpsf,thresh, mkiraf=True):
    """
    Performs sextracting by daofind and photometry by daophot 
    filename is WITHOUT suffix .fits
    fwhmpsf is rough estimation of typical frame FWHM
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
    iraf.dNb_cuts = Nb_cutsigiphot
    iraf.daophot

    iraf.unlearn('phot')
    iraf.unlearn('datapars')
    iraf.unlearn('photpars')
    iraf.unlearn('daopars')

    iraf.datapars.epadu = 1
    iraf.datapars.readnoi = 6
    iraf.datapars.datamin = "INDEF"
    iraf.datapars.datamax = "50000"

    # fwhm tarot 1.5   oaj 3.5
    iraf.datapars.fwhm = fwhmpsf

    for filename in filelist:

        IMstatres = iraf.imstat(filename,Stdout=1)
        IMmean = IMstatres[1].split()[2]
        iraf.datapars.sigma = math.sqrt(float(IMmean) * float(iraf.datapars.epadu) + float(iraf.datapars.readnoi)**2) / float(iraf.datapars.epadu)

        print("--- performing daophot sextracting ---")
        iraf.daofind(filename,output="default",verify="no", verbose="no", threshold=thresh)
        iraf.datapars.datamax = "INDEF"
        print("--- performing daophot photometry ---")
        iraf.daophot.phot(image=filename,coords="default", output="default", interactive="no", sigma="INDEF", airmass="AIRMASS", exposure="EXPOSURE", filter="FILTER", obstime="JD", calgorithm="gauss", verify="no", verbose="no")
    

def select_good_stars(filelist,limiting_mag_err):
# Performs selection of stars without INDEF in magnitude or magnitude error and with the flag "NoError" inside *.mag.1 file
# Saves into *.magfiltered file in x-y coordinates
# filename is WITHOUT suffix .fits

    for filename in filelist:
        path, filename2 = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''
        
        magfile = filename2 + '.mag.1'
        resmaggile = filename2 + ".magfiltered"
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
                print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
            else:
                #print(merr)
                if (float(merr) < limiting_mag_err):
                    #print(xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                    f2.write(xpos+" "+ypos+" "+mag+" "+merr+"\n")                
                else:
                    print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)

            lajna = f1.readline()

        f1.close()
        f2.close()


def convert_xy_radec(filelist, soft='sextractor'):
    """
    Performs pyraf transformation of x-y into RA-DEC coordinates
    filename is WITHOUT suffix .fits
    Input is the *.magfiltered file from select_good_stars() function
    Output is the *.magwcs file
    """
    for filename in filelist:
        path, filename_ext = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        # Get rid of the extension to keep only the name
        filename2 = filename_ext.split('.')[0]

        magfilewcs = folder+filename2 + ".magwcs"

        if soft == 'iraf':
            from pyraf import iraf

            magfile = folder+filename2 + ".magfiltered"
            iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")
            data1 = ascii.read(magfile, names=['Xpos', 'Ypos', 'Mag_aper', 'Mag_err_aper' ])
            data2 = ascii.read(magfilewcs, names=['RA', 'DEC', 'Mag_aper', 'Mag_err_aper' ])
            data = Table([data1['Xpos'],data1['Ypos'], data2['RA'], data2['DEC'], data2['Mag_aper'], data2['Mag_err_aper'], [filename]*len(data1)], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_inst', 'Magerr_inst', 'filenames'])
 
        elif soft == 'sextractor':
            sources = ascii.read(folder + filename2 + '_SourcesDet.cat', format='sextractor')
            header = fits.getheader(filename)
            w = wcs.WCS(header)
            ra, dec = w.wcs_pix2world(sources['X_IMAGE'], sources['Y_IMAGE'], 1)
            filenames = [filename] * len(ra)
            data = Table([sources['X_IMAGE'],  sources['Y_IMAGE'], ra,dec, sources['MAG_AUTO'], sources['MAGERR_AUTO'], filenames], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_isnt', 'Magerr_inst', 'filenames'])

        data.write(magfilewcs, format='ascii.commented_header', overwrite=True)
        #data4=data['RA', 'DEC']
        #data4.write(magfilewcs+'2',format='ascii.commented_header', overwrite=True)


def crosscheck_with_catalogues(image_table, radius, catalogs=['I/284/out', 'I/345/gaia2', 'II/349/ps1', 'I/271/out'], Nb_cuts=(1,1)):
#def crosscheck_with_catalogues(image_table, radius, catalogs=['I/284/out'], Nb_cuts=(1,1)):
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
    transients_out_list = []
    for i, filename in enumerate(image_table['filenames']):
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
        else:
            split_file = filename2.split('_')
            if len(split_file[-1]) == 2:
                idx = 3
            elif len(split_file[-1]) == 3:
                idx = 4
            original_name = folder + filename2[:-idx]

        transients_out_list.append(original_name + ".oc")

        header = fits.getheader(filename)

        # Get pixel scale in degrees
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

        # Load detected sources in astropy table
        detected_sources = ascii.read(magfilewcs, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames'])
        detected_sources['quadrant'] = [image_table['quadrant'][i]]*len(detected_sources)

        # Transform X_pos and Y_pos to original image in case it was split
        header2 = fits.getheader(original_name + extension)
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
        detected_sources['Xpos_quad'] = detected_sources['Xpos'] 
        detected_sources['Ypos'] = detected_sources['Ypos'] + Naxis11 * int(index_i)
        detected_sources['Ypos_quad'] = detected_sources['Ypos'] 

        if i == 0:
            detected_sources_tot = deepcopy(detected_sources)
        else:
            detected_sources_tot = vstack([detected_sources_tot, detected_sources])
    # Add units
    detected_sources_tot['_RAJ2000'] *= u.deg
    detected_sources_tot['_DEJ2000'] *= u.deg
    # Add index for each source
    detected_sources_tot['idx'] = np.arange(len(detected_sources_tot))
    
    # Initialise candidates with all detected sources
    candidates = deepcopy(detected_sources_tot)
    candidates.write('test0.dat', format='ascii.commented_header', overwrite=True)
    print ('\nCrossmatching sources with catalogs.')
    print ('Radius used for crossmatching with catalogs: %.2f arcseconds\n' % (radius*pixScale*3600))
    for catalog in catalogs:
        print (catalog, len(candidates))
        # Use Xmatch to crossmatch with catalog
        crossmatch = run_xmatch(candidates[:10], catalog, radius*pixScale*3600)

        crossmatch.write('test.dat', format='ascii.commented_header', overwrite=True)
        # Initialise flag array. 0: unknown sources / 1: known sources
        flag = np.zeros(len(candidates))
        # Do not consider duplicates
        referenced_star_idx = np.unique(crossmatch['idx'])
        #print (referenced_star_idx)
        #print (candidates)
        #print (np.max(referenced_star_idx))
        #ref_idx = 
        # Set flag indexes to 1 for detected sources associated to a star
        flag[np.array(referenced_star_idx)] = 1
        #print (len(candidates)) 
        #print (flag)
        # Table for candidates
         
        candidates = candidates[flag == 0]
        #print (len(candidates))
        # Update indexes
        candidates['idx'] = np.arange(len(candidates))
        print ('%d/%d candidates left after crossmatching with %s' % (len(candidates),len(detected_sources), cat_dict[catalog]))
        if (len(candidates) == 0):
            break
        
    transients = np.unique(transients_out_list)[0]
    # Write candidates file
    candidates.write(transients, format='ascii.commented_header', overwrite=True)
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
                        type=float,
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
                        choices=['TRE','TCA','TCH','OAJ','Lisniky-AZT8','UBAI-T60S','UBAI-T60N'],
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




    args = parser.parse_args()

    Nb_cuts = (args.quadrants,args.quadrants)
    
    # Load config files for a given telescope
    config = load_config(args.telescope)

    image_table = astrometric_calib(args.filename, config, Nb_cuts=Nb_cuts,soft=args.doAstrometry, verbose=args.verbose)

    if args.soft == 'iraf':
        clean_folder(image_table['filenames'])
        get_photometry(image_table['filenames'], args.FWHM, args.threshold)
        select_good_stars(image_table['filenames'], args.mag_err_cut)
    elif args.soft == 'sextractor':
        sextractor(image_table['filenames'], args.FWHM, args.threshold, args.telescope,config,verbose=args.verbose)

    convert_xy_radec(image_table['filenames'],soft=args.soft)
    total_candidates = crosscheck_with_catalogues(image_table,args.radius_crossmatch, Nb_cuts=Nb_cuts)
    #check_moving_objects(args.filename, total_candidates)
   
    #total_candidates = ascii.read('total_candidates.dat', names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames', 'idx' ,'quadrant'])
    #total_candidates_calib = phot_calib(total_candidates, args.telescope, radius=args.radius_crossmatch,doPlot=True)
    #total_candidates_calib = ascii.read('tot_cand2.dat', names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames', 'idx' ,'quadrant', 'mag_calib', 'mag_calib_err', 'magsys', 'filter_cat', 'filter_DB'])

    # If both arguments VOE_path and owncloud_path are provided
    # Send candidates to database
    if args.VOE_path and args.owncloud_path:
        send_data2DB(args.filename, total_candidates_calib, Nb_cuts, args.owncloud_path, args.VOE_path, "utilsDB/usrpwd.json",debug=True)



