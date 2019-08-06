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
#   python gmadet.py --filename /folder/image.fits --FWHM 1.5 --threshold 4 --radius_crossmatch 3 --soft iraf
#   python gmadet.py --filename /folder/image.fits --FWHM 3.5 --threshold 4 --radius_crossmatch 3 --soft sextractor
# #   python gmadet.py --filename /folder/image.fits --FWHM 4 --threshold 10 --radius_crossmatch 3 --soft sextractor

import sys, subprocess, glob, math, shutil, os
import argparse

from catalogues import *
from astromatic import psfex, scamp

from astropy.io import ascii, fits
from astropy.table import vstack, Table

from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
import voeventparse as vp
import json

def clean_folder(filename):
    """ Remove output files from previous iraf run. No need for sextractor  """

    types = ('*coo.*', '*mag.*', '*.magwcs', '*.magfiltered*')
    files2delete = []
   
    path = os.path.split(filename)
    if path[0]:
        folder = path[0] + '/'
    else:
        folder = ''

    for f in types:
        files2delete.extend(glob.glob(folder+f))
    for f in files2delete:
        os.remove(f)

def astrometric_calib(filename):
    """perform astrometric calibration"""
    scamp(filename)

def psfex(filename, telescope):
    """Run psfex to estimate PSF FWHM"""

    subprocess.call(['sex', '-c', 'config/%s/sourcesdet.sex' % telescope, filename, '-SEEING_FWHM', str(fwhmpsf), '-DETECT_THRESH', str(thresh), '-PARAMETERS_NAME', 'config/%s/sourcesdet.param' % telescope])


def sextractor(filename, fwhmpsf, thresh, telescope):
    """Run sextractor """

    subprocess.call(['sex', '-c', 'config/%s/sourcesdet.sex' % telescope, filename, '-SEEING_FWHM', str(fwhmpsf), '-DETECT_THRESH', str(thresh), '-PARAMETERS_NAME', 'config/%s/sourcesdet.param' % telescope])



def get_photometry(filename,fwhmpsf,thresh, mkiraf=True):
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


def select_good_stars(filename,limiting_mag_err):
# Performs selection of stars without INDEF in magnitude or magnitude error and with the flag "NoError" inside *.mag.1 file
# Saves into *.magfiltered file in x-y coordinates
# filename is WITHOUT suffix .fits
    list_of_magfiles = glob.glob("./*.mag.*")
    magfile = list_of_magfiles[0]
    resmaggile = filename+".magfiltered"
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


def convert_xy_radec(filename, soft='sextractor'):
    """
    Performs pyraf transformation of x-y into RA-DEC coordinates
    filename is WITHOUT suffix .fits
    Input is the *.magfiltered file from select_good_stars() function
    Output is the *.magwcs file
    """
    magfilewcs = filename+".magwcs"

    if soft == 'iraf':
        from pyraf import iraf

        magfile = filename+".magfiltered"
        iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")
        #shutil.remove(magfile)

    elif soft == 'sextractor':
        sources = ascii.read('sourcesdet.cat', format='sextractor')
        header = fits.getheader(filename)
        w = wcs.WCS(header)
        ra, dec = w.wcs_pix2world(sources['X_IMAGE'], sources['Y_IMAGE'], 0)
        data = Table([sources['X_IMAGE'],  sources['Y_IMAGE'], ra,dec, sources['MAG_APER'], sources['MAGERR_APER']], names=['Xpos', 'Ypos', 'RA', 'DEC', 'Mag_aper', 'Mag_err_aper'])
        data.write(magfilewcs, format='ascii.commented_header')

def crosscheck_with_catalogues(filename, radius, catalogs=['I/284/out', 'I/345/gaia2', 'II/349/ps1', 'I/271/out']):
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

    magfilewcs = filename+".magwcs"
    transients = filename+".oc"

    header = fits.getheader(filename)

    # Get pixel size
    #pixSize = abs(header['CDELT1'])
    pixSize = abs(header['CD1_1'])
    print (pixSize)

    # Load detected sources in astropy table
    #detected_sources = ascii.read(magfilewcs, names=['_RAJ2000','_DEJ2000', 'mag_int', 'mag_inst_err'])
    detected_sources = ascii.read(magfilewcs, names=['Xpos','Ypos','_RAJ2000','_DEJ2000', 'mag_inst', 'mag_inst_err'])
    # Add units
    detected_sources['_RAJ2000'] *= u.deg
    detected_sources['_DEJ2000'] *= u.deg
    # Add index for each source
    detected_sources['idx'] = np.arange(len(detected_sources))
    
    # Initialise candidates with all detected sources
    candidates = detected_sources

    for catalog in catalogs:
        print (catalog)
        # Use Xmatch to crossmatch with catalog
        crossmatch = xmatch(candidates, catalog, radius*pixSize*3600)
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

    # Write candidates file
    candidates.write(transients, format='ascii.commented_header', overwrite=True)

    return candidates

def check_moving_objects(filename, candidates):
    """
    Crossmatch the list of candidates with moving objects using SkyBoT
    
    """
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

    return candidates_out


def create_subimages(filename, candidates, owncloud_path, VOE_path, usrpwd_path, coords_type='pix', corner_cut=32):
    """Create sub-images centered on OT position"""

    from create_OT_subimage import make_sub_image

    # Load data and header
    data, header = fits.getdata(filename, header = True)
    dateObs = header['DATE-OBS']
    Tstart = Time(dateObs, format='fits', scale='utc')
    Tend = Tstart + TimeDelta(float(header['EXPOSURE']), format='sec')
    Airmass = header['AIRMASS']

    # Do not consider candidates found in the image edge
    imsize = data.shape
    mask = (candidates['Xpos'] > corner_cut) & \
            (candidates['Ypos'] > corner_cut) & \
            (candidates['Xpos'] < imsize[1] - corner_cut) & \
            (candidates['Ypos'] < imsize[0] - corner_cut)
    candidates_cut = candidates[mask]

    # Get information about the current alert from the xml file containing observation plan
    with open(VOE_path,'rb') as f:
        obsplan = vp.load(f)

    dict_event={}
    dict_event['event_type'] = obsplan.find(".//Param[@name='Event_type']").attrib['value']
    dict_event['event_name'] = obsplan.find(".//Param[@name='Event_ID']").attrib['value']
    dict_event['event_status'] = obsplan.find(".//Param[@name='Event_status']").attrib['value']
    dict_event['revision'] = obsplan.find(".//Param[@name='Revision']").attrib['value']
    dict_event['telescope'] = obsplan.find(".//Param[@name='Name_tel']").attrib['value']

    # Get user email adress and password to login in https://grandma-fa-interface.lal.in2p3.fr
    with open(usrpwd_path) as f:
        usrpwd = json.load(f)


    # Set up the output repository path to store sub-images
    outputDir = owncloud_path + '/' + dict_event['event_type'] + '/' + \
            dict_event['event_name'] + '/' + dict_event['event_status'] + \
            '_' + dict_event['revision'] + '/OTs/'

    # Create a sub image centered on each candidate found, and gather information
    tile_id_list = [1] * len(candidates_cut)
    filter_list = ['Clear'] * len(candidates_cut)
    Tstart_list = [Tstart.fits] * len(candidates_cut)
    Tend_list = [Tend.fits] * len(candidates_cut)
    Airmass_list = [Airmass] * len(candidates_cut)
    Fits_path = []
    for i, row in enumerate(candidates_cut):
        name = dict_event['telescope'] + '_' + \
                str(round(float(row['_RAJ2000']),5)) + '_' + \
                str(round(float(row['_DEJ2000']),5)) + '_' + \
                dateObs + '.fits.gz'
        #name = 'test' + str(i) + '.fits.gz'
        Fits_path.append(name)
        if coords_type == 'world':
            OT_coords = [row['_RAJ2000'], row['_DEJ2000']]
        elif coords_type == 'pix':
            OT_coords = [row['Ypos'], row['Xpos']]
        make_sub_image(filename, OT_coords, coords_type=coords_type,
                       output_name=outputDir+name, size=[128,128])

    alias = ['new'] * len(candidates_cut)
    new = [1] * len(candidates_cut)
    tile_id_list = [2] * len(candidates_cut)
    RA_list = candidates_cut['_RAJ2000']
    Dec_list = candidates_cut['_DEJ2000']
    filter_list = ['Clear'] * len(candidates_cut)
    Tstart_list = [Tstart.fits] * len(candidates_cut)
    Tend_list = [Tend.fits] * len(candidates_cut)
    Mag_list = candidates_cut['mag_inst']+20
    Mag_err_list = candidates_cut['mag_inst_err']
    Magsys_list = ['Instrumental'] * len(candidates_cut)
    Airmass_list = [Airmass] * len(candidates_cut)

    candidates_2DB = Table([alias,new,tile_id_list,RA_list,Dec_list,filter_list,Tstart_list,Tend_list,Mag_list,Mag_err_list,Magsys_list,Airmass_list,Fits_path], names=['alias','new','tile_id','RA','DEC','filter','Tstart','Tend','Magnitude','Magnitude_error','Magsys','Airmass','fits_name'])

    
    # Set url to report tile or galaxy observations
    url = "https://grandma-fa-interface.lal.in2p3.fr/obs_report_OT.php"
    #url = "http://localhost/test2/obs_report_OT.php"

    # Loop over the observations
    for i in range(len(candidates_2DB)):
        data2DB={}
        for col in candidates_2DB.colnames:
            data2DB[col] = candidates_2DB[col][i]

        # Add obsplan info to data dictionary
        for key, value in dict_event.items():
            data2DB[key] = value

        # Add username and password to data dictionary
        for key, value in usrpwd.items():
            data2DB[key] = value

        # Add compulsory keys
        data2DB["method"] = "POST"
        data2DB["submit"] = "ok"

        response = requests.post(url, data=data2DB)

        print ('\nDEBUG:\n')
        print ('Data sent to DB:')
        print (data2DB)
        print ('\n\n')
        print ('Request response text:')
        print (response.text)
        print ('\n\n')
        print ('Request response status code:')
        print (response.status_code)
        print ('\n\n')
        print ('Request response history:')
        print (response.history)


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
                        required=True,
                        type=str,
                        help='Local path to the owncloud')

    parser.add_argument('--VOE_path',
                        dest='VOE_path',
                        required=True,
                        type=str,
                        help='Path + filename of the VoEvent containing the observation plan')

    parser.add_argument('--usrpwd_path',
                        dest='usrpwd_path',
                        required=True,
                        type=str,
                        help='Path + filename of the json file containing authentification information for https://grandma-fa-interface.lal.in2p3.fr')

    args = parser.parse_args()


    #astrometric_calib(args.filename)
    #psfex(args.filename)
    if args.soft == 'iraf':
        clean_folder(args.filename)
        get_photometry(args.filename, args.FWHM, args.threshold)
        select_good_stars(args.filename, args.mag_err_cut)
    elif args.soft == 'sextractor':
        sextractor(args.filename, args.FWHM, args.threshold, args.telescope)

    convert_xy_radec(args.filename,soft=args.soft)
    print ('Crossmatch with catalogs')
    candidates = crosscheck_with_catalogues(args.filename,args.radius_crossmatch)
    #check_moving_objects(args.filename, candidates)

    create_subimages(args.filename, candidates, args.owncloud_path, args.VOE_path, args.usrpwd_path)



