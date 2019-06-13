#! /usr/bin/env python

# Python-Pyraf module for GRANDMA detection of Optical candidates
# Authors: Martin Blazek, Granada, Spain, alf@iaa.es
#          David Corre, Orsay, France, corre@lal.in2p3.fr
# v1.1, last modified 2019 June
# Good luck everyone with pyraf installation
# Input arguments are filename without fits suffix, typical fwhm, sextracting threshold and maximal distance for catalogue crosschecking in degrees
# 
# Example:
#   python3 baf.py tarot 1.5 4 0.000907203
#   python3 baf.py oaj 3.5 4 0.0001543390

import sys
import math
import glob
from catalogues import *

from pyraf import iraf
from pyraf.iraf import daophot

# Parameter to get rid of the faint stars with high magnitude error
magnitude_error_threshold = 0.5

# Parameter for catalogue crosschecking, maximal allowed distance between sextracted and catalogue position in degrees
# Tarot 1px = 0.000907203 deg
# OAJ 1px = 1.543390792967E-04 deg
# allowed_crosscheck_radius = 3*0.000907203       # Tarot
# allowed_crosscheck_radius = 3*0.0001543390      # OAJ

allowed_crosscheck_radius = sys.argv[4]

def get_photometry(filename,fwhmpsf,THRESH):
# Performs sextracting by daofind and photometry by daophot 
# filename is WITHOUT suffix .fits
# fwhmpsf is rough estimation of typical frame FWHM
# THRESH is signal-to-noise ratio, typically 4 or so
# Outputs are two files *.coo.1 and *.mag.1 with x-y coordinates

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
    iraf.daofind(filename,output="default",verify="no", verbose="no", threshold=THRESH)
    iraf.datapars.datamax = "INDEF"
    print("--- performing daophot photometry ---")
    iraf.daophot.phot(image=filename,coords="default", output="default", interactive="no", sigma="INDEF", airmass="AIRMASS", exposure="EXPOSURE", filter="FILTER", obstime="JD", calgorithm="gauss", verify="no", verbose="no")
    


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
            print(merr)
            if (float(merr) < limiting_mag_err):
                print(xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)
                f2.write(xpos+" "+ypos+" "+mag+" "+merr+"\n")                
            else:
                print("XXXXXX "+xpos+" "+ypos+" "+mag+" "+merr+" "+errmessage)

        lajna = f1.readline()

    f1.close()
    f2.close()

def convert_xy_radec(filename):
# Performs pyraf transformation of x-y into RA-DEC coordinates
# filename is WITHOUT suffix .fits
# Input is the *.magfiltered file from select_good_stars() function
# Output is the *.magwcs file
    magfile = filename+".magfiltered"
    magfilewcs = filename+".magwcs"
    iraf.wcsctran(input=magfile, output=magfilewcs, image=filename, inwcs="physical", outwcs="world")


def crosscheck_with_catalogues(filename,degrad):
# Performs crosscheck with USNO B1.0 catalogue with *.magwcs
# filename is WITHOUT suffix .fits and maximal allowed difference degrad is in degrees
# Input file is *.magwcs and the output is the list of the stars *.oc which were not identified in the catalogue
    magfilewcs = filename+".magwcs"
    transients = filename+".oc"
    
    f3 = open(magfilewcs,"r")
    lajna = f3.readline()
    lajna = f3.readline()
    lajna = f3.readline()
    lajna = f3.readline()
    
    f4 = open(transients,"w")
    counter = 1
    while lajna:
        #print(lajna)
        ra = lajna.split()[0]
        dec = lajna.split()[1]
        mag = lajna.split()[2]
        merr = lajna.split()[3]
        try:
            wohoo = USNO_B1_query(ra,dec,float(eval(degrad)))
            #print(wohoo)
        except:
            print("New transient "+str(counter)+" at "+ra+" "+dec)
            counter=counter+1
            f4.write(ra+" "+dec+" "+mag+" "+merr+"\n")
        lajna=f3.readline()
    
    f3.close()
    f4.close()
    


get_photometry(sys.argv[1],sys.argv[2],sys.argv[3])
select_good_stars(sys.argv[1],magnitude_error_threshold)
convert_xy_radec(sys.argv[1])
crosscheck_with_catalogues(sys.argv[1],allowed_crosscheck_radius)






