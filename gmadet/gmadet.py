#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Python-Pyraf module for GRANDMA detection of Optical candidates
# Authors: Martin Blazek, Granada, Spain, alf@iaa.es
#          David Corre, Orsay, France, corre@lal.in2p3.fr
# v1.1, last modified 2019 June
# Good luck everyone with pyraf installation
# Input arguments are filename without fits suffix, typical fwhm, sextracting threshold and maximal distance for catalogue crosschecking in degrees
# 
# Example:
#   python gmadet.py --filename /folder/image.fits --FWHM 1.5 --threshold 4 --radius_crossmatch 3 --soft iraf --telescope TRE --owncloud_path /home/corre/ownCloud --VOE_path /home/corre/ownCloud/GW/Test_Reporting_Obs/Initial_0/Obs_request/GRANDMA20190505_GWTest_Reporting_Obs_ShAO-T60_a.xml --quadrants 4 --doAstrometry True
#   python gmadet.py --filename /folder/image.fits --FWHM 3.5 --threshold 4 --radius_crossmatch 3 --soft sextractor --telescope TRE --owncloud_path /home/corre/ownCloud --VOE_path /home/corre/ownCloud/GW/Test_Reporting_Obs/Initial_0/Obs_request/GRANDMA20190505_GWTest_Reporting_Obs_ShAO-T60_a.xml --quadrants 1 --doAstrometry False

import sys, subprocess, glob, math, shutil, os
import argparse
import warnings

from catalogues import *
from phot_calibration import phot_calib
from utils import send_data2DB
from astromatic import psfex, scamp

from astropy.io import ascii, fits
from astropy.table import vstack, Table, Column

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

from copy import deepcopy

warnings.simplefilter(action='ignore', category=FutureWarning)

def clean_folder(filelist):
    """ Remove output files from previous iraf run. No need for sextractor  """

    types = ('*coo.*', '*mag.*', '*.magwcs', '*.magfiltered*')
    files2delete = []
    for filename in  filelist:
        path = os.path.split(filename)
        if path[0]:
            folder = path[0] + '/'
        else:
            folder = ''

        for f in types:
            files2delete.extend(glob.glob(folder+f))
            files2delete.extend(glob.glob(f))

    files2delete = np.unique(files2delete)
    for f in files2delete:
        os.remove(f)

def astrometric_calib(filename, Nb_cuts=(2,2), perform=False):
    """perform astrometric calibration"""
    from astro_calibration import cut_image, perform_astrometry

    path, filename2 = os.path.split(filename)
    if path:
        folder = path + '/'
    else:
        folder = ''
    cut_image(filename, Nb_cuts = Nb_cuts)

    filelist = []
    quadrant = []
    
    if perform:
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

        scaleLow = 0.7 * pixScale * 3600
        scaleHigh = 1.3 * pixScale * 3600
        radius = max(header['NAXIS1']*pixScale, header['NAXIS2']*pixScale)
        if Nb_cuts[0] > 1 or Nb_cuts[1] > 1:
            index = 0
            for i in range(Nb_cuts[0]):
                for j in range(Nb_cuts[1]):
                    index+=1
                    filein = folder + "D%d_" % (index) + filename2
                    filelist.append(filein)
                    perform_astrometry(filein,radius=radius,scaleLow=scaleLow,scaleHigh=scaleHigh)
                    quadrant.append('Q%d_%d_%d' % (index,i,j))
        else:
            filelist.append(folder + "D1_" + filename2)
            perform_astrometry(folder + "D1_" + filename2,radius=radius,scaleLow=scaleLow,scaleHigh=scaleHigh)
            quadrant.append('Q1_0_0')
    else:
        filelist.append(folder + "D1_" + filename2)
        quadrant.append('Q1_0_0')

    image_table = Table([filelist, quadrant], names=['filenames', 'quadrant']) 
    return image_table

def psfex(filename, telescope):
    """Run psfex to estimate PSF FWHM"""

    subprocess.call(['sex', '-c', 'config/%s/sourcesdet.sex' % telescope, filename, '-SEEING_FWHM', str(fwhmpsf), '-DETECT_THRESH', str(thresh), '-PARAMETERS_NAME', 'config/%s/sourcesdet.param' % telescope])


def sextractor(filelist, fwhmpsf, thresh, telescope):
    """Run sextractor """

    for filename in filelist:
        path, filename2 = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''
        print (filename2)
        print (filename2.split('.fits')[0])
        subprocess.call(['sex', '-c', 'config/%s/sourcesdet.sex' % telescope, filename, '-SEEING_FWHM', str(fwhmpsf), '-DETECT_THRESH', str(thresh), '-PARAMETERS_NAME', 'config/%s/sourcesdet.param' % telescope, '-CATALOG_NAME', folder + 'sourcesdet_%s.cat' % (filename2.split('.fits')[0])])



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
    iraf.digiphot
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
    


def run_psfex(filename):
    psfex(filename) 
    psffilename = filename.split('.')[0] + '.psf.fits'


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
        path, filename2 = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        magfilewcs = filename + ".magwcs"

        if soft == 'iraf':
            from pyraf import iraf

            magfile = filename2 + ".magfiltered"
            iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")
            data1 = ascii.read(magfile, names=['Xpos', 'Ypos', 'Mag_aper', 'Mag_err_aper' ])
            data2 = ascii.read(magfilewcs, names=['RA', 'DEC', 'Mag_aper', 'Mag_err_aper' ])
            data = Table([data1['Xpos'],data1['Ypos'], data2['RA'], data2['DEC'], data2['Mag_aper'], data2['Mag_err_aper'], [filename]*len(data1)], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_inst', 'Magerr_inst', 'filenames'])
 
        elif soft == 'sextractor':
            sources = ascii.read(folder + 'sourcesdet_%s.cat' % (filename2.split('.fits')[0]), format='sextractor')
            header = fits.getheader(filename)
            w = wcs.WCS(header)
            ra, dec = w.wcs_pix2world(sources['X_IMAGE'], sources['Y_IMAGE'], 1)
            filenames = [filename] * len(ra)
            data = Table([sources['X_IMAGE'],  sources['Y_IMAGE'], ra,dec, sources['MAG_AUTO'], sources['MAGERR_AUTO'], filenames], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_isnt', 'Magerr_inst', 'filenames'])

        data.write(magfilewcs, format='ascii.commented_header', overwrite=True)
        data4=data['RA', 'DEC']
        data4.write(magfilewcs+'2',format='ascii.commented_header', overwrite=True)


def crosscheck_with_catalogues(image_table, radius, catalogs=['I/284/out', 'I/345/gaia2', 'II/349/ps1', 'I/271/out'], Nb_cuts=(1,1)):
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
        path, filename2 = os.path.split(filename)
        if path:
            folder = path + '/'
        else:
            folder = ''

        magfilewcs = filename+".magwcs"

        split_file = filename2.split('_')
        if len(split_file[0]) == 2:
            idx = 3
        elif len(split_file[0]) == 3:
            idx = 4
        transients_out_list.append(folder + filename2[idx:] + ".oc")

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
        #detected_sources = ascii.read(magfilewcs, names=['_RAJ2000','_DEJ2000', 'mag_int', 'mag_inst_err'])
        detected_sources = ascii.read(magfilewcs, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames'])
        detected_sources['quadrant'] = [image_table['quadrant'][i]]*len(detected_sources)

        # Transform X_pos and Y_pos to original image in case it was split
        header2 = fits.getheader(folder + filename2[idx:])
        Naxis1 = float(header2['NAXIS1'])
        Naxis2 = float(header2['NAXIS2'])
        Naxis11 = int(Naxis1/Nb_cuts[0])
        Naxis22 = int(Naxis2/Nb_cuts[1])
        quad, index_i, index_j = image_table['quadrant'][i].split('_')
        quad = quad[1:]
        detected_sources['Xpos'] = detected_sources['Xpos'] + Naxis22 * int(index_j)
        detected_sources['Ypos'] = detected_sources['Ypos'] + Naxis11 * int(index_i)

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
    candidates = detected_sources_tot
    print ('\nCrossmatching sources with catalogs.')
    print ('Radius used for crossmatching with catalogs: %.2f arcseconds\n' % (radius*pixScale*3600))
    for catalog in catalogs:
        #print (catalog)
        # Use Xmatch to crossmatch with catalog
        crossmatch = xmatch(candidates, catalog, radius*pixScale*3600)
        # Initialise flag array. 0: unknown sources / 1: known sources
        flag = np.zeros(len(candidates))
        # Do not consider duplicates
        referenced_star_idx = np.unique(crossmatch['idx'])
        # Set flag indexes to 1 for detected sources associated to a star
        flag[referenced_star_idx] = 1
        
        # Table for candidates
        candidates = candidates[flag == 0]
        # Update indexes
        candidates['idx'] = np.arange(len(candidates))
        print ('%d/%d candidates left after crossmatching with %s' % (len(candidates),len(detected_sources), cat_dict[catalog]))
        if (len(candidates) == 0):
            break
        
    transients = np.unique(transients_out_list)[0]
    # Write candidates file
    candidates.write(transients, format='ascii.commented_header', overwrite=True)
    """
        # If some candidates found in one quadrant, add them to the total 
        if len(candidates) > 0:
            counter += 1
            candidates2 = candidates
            quadrantcol = Column(np.array([image_table['quadrant'][i]]*len(candidates)), name='quadrant')
            candidates2.add_column(quadrantcol)
            if counter == 1:
                total_candidates = deepcopy(candidates2)
            elif counter > 1:
                total_candidates = vstack([total_candidates, candidates2])

    # Write all the candidates found in each subimage
    total_candidates.write('total_candidates.dat', format='ascii.commented_header', overwrite=True)
    """
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
                        required=True,
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
                        help='Path + filename of the VoEvent containing the observation plan')

    parser.add_argument('--quadrants',
                        dest='quadrants',
                        required=True,
                        type=int,
                        help='Number of quadrants the image is divided')

    parser.add_argument('--doAstrometry',
                        dest='doAstrometry',
                        required=True,
                        choices=['True', 'False'],
                        type=str,
                        help='Whether to perform astrometric calibration')
    args = parser.parse_args()
 
    Nb_cuts = (args.quadrants,args.quadrants)
    
    if args.doAstrometry == 'True':
        doAstrometry = True
    elif args.doAstrometry:
        doAstrometry = False

    image_table = astrometric_calib(args.filename, Nb_cuts=Nb_cuts,perform=doAstrometry)
    #psfex(args.filename)
    if args.soft == 'iraf':
        clean_folder(image_table['filenames'])
        get_photometry(image_table['filenames'], args.FWHM, args.threshold)
        select_good_stars(image_table['filenames'], args.mag_err_cut)
    elif args.soft == 'sextractor':
        sextractor(image_table['filenames'], args.FWHM, args.threshold, args.telescope)

    convert_xy_radec(image_table['filenames'],soft=args.soft)
    total_candidates = crosscheck_with_catalogues(image_table,args.radius_crossmatch, Nb_cuts=Nb_cuts)
    #check_moving_objects(args.filename, total_candidates)
   
    #total_candidates = ascii.read('total_candidates.dat', names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames', 'idx' ,'quadrant'])
    total_candidates_calib = phot_calib(total_candidates, args.telescope, radius=args.radius_crossmatch,doPlot=False)
    #total_candidates_calib = ascii.read('tot_cand2.dat', names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err', 'filenames', 'idx' ,'quadrant', 'mag_calib', 'mag_calib_err', 'magsys', 'filter_cat', 'filter_DB'])
    send_data2DB(args.filename, total_candidates_calib, Nb_cuts, args.owncloud_path, args.VOE_path, "utilsDB/usrpwd.json",debug=True)



