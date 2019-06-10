# -*- coding: utf-8 -*-
  
"""Cross match candidates with catalogs."""

from astropy.io import ascii
from astropy.table import Table

import sys, re, os
import numpy as np
import pylab
import json
import requests
import h5py

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib   


def gaia_query(ra_deg, dec_deg, rad_deg, maxmag=20,
               maxsources=-1):
    """
    Query Gaia DR2 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RAJ2000', 'e_RAJ2000',
                             'DEJ2000','e_DEJ2000','Dup','Gmag','e_Gmag', 'o_Gmag',
                             'BPmag', 'e_BPmag', 'o_BPmag', 'RPmag', 'e_RPmag', 'o_RPmag'],
    #                column_filters={"gmag":
    #                                ("<%f" % maxmag),
    #                               "imag":
    #                                ("<%f" % maxmag)},
                    row_limit = maxsources)
    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')

    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="I/345/gaia2")[0]

def sdss_query(ra_deg, dec_deg, rad_deg, maxmag=20,
               maxsources=-1):
    """
    Query SDSS DR12 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RA_ICRS', 'e_RA_ICRS',
                             'DE_ICRS','e_DE_ICRS','umag','e_umag','gmag', 'e_gmag',
                             'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag',
                             'zsp','spCl', 'subCl'],
    #                column_filters={"gmag":
    #                                ("<%f" % maxmag),
    #                               "imag":
    #                                ("<%f" % maxmag)},
                    row_limit = maxsources)
    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="V/147/sdss12")[0]

def _2MASS_query(ra_deg, dec_deg, rad_deg, maxmag=20,
               maxsources=-1):
    """
    Query 2MASS @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['2MASS', 'RAJ2000', 'DEJ2000',
                             'Jmag','e_Jmag','Jsnr',
                             'Hmag','e_Hmag','Hsnr',
                             'Kmag','e_Kmag','Ksnr',
                             'dup', 'Ndet'],
    #                column_filters={"gmag":
    #                                ("<%f" % maxmag)
    #                               "imag":
    #                                ("<%f" % maxmag)},
                    row_limit = maxsources)
    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/246")[0]

def USNO_B1_query(ra_deg, dec_deg, rad_deg, maxmag=20,
               maxsources=-1):
    """
    Query USNO_B1 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['USNO-B1.0', 'RAJ2000', 'e_RAJ2000','DEJ2000', 'e_DEJ2000',
                             'B1mag','R1mag','B2mag',
                             'R2mag','Imag', 'Ndet'],
    #                column_filters={"gmag":
    #                                ("<%f" % maxmag)
    #                               "imag":
    #                                ("<%f" % maxmag)},
                    row_limit = maxsources)
    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="I/284/out")[0]

def USNO_A2_query(ra_deg, dec_deg, rad_deg, maxmag=20,
               maxsources=-1):
    """
    Query USNO_A2 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['USNO-A2.0', 'RAJ2000', 'DEJ2000',
                             'Bmag','Rmag'],
    #                column_filters={"gmag":
    #                                ("<%f" % maxmag)
    #                               "imag":
    #                                ("<%f" % maxmag)},
                    row_limit = maxsources)
    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="I/252/out")[0]

def get_glade():
    if not os.path.isdir('catalogs/'):
        os.makedirs('catalogs/')

    catalogFile = os.path.join('catalogs/', "glade.hdf5")
        
    if not os.path.isfile(catalogFile):
        # Unset row limits when querying Vizier
        Vizier.ROW_LIMIT = -1
        cat, = Vizier.get_catalogs('VII/281/glade2')
        ra, dec = cat["RAJ2000"], cat["DEJ2000"]
        distmpc, z, flag1 = cat["Dist"], cat["z"], cat["Flag1"]
        magb, BMAG = cat["Bmag"], cat["BMAG"]
        Jmag, Hmag, Kmag = cat["Jmag"], cat["Hmag"], cat['Kmag']
        flag2, flag3 =  cat["Flag2"],  cat["Flag3"]
        # Keep track of galaxy identifier
        GWGC, PGC, HyperLEDA = cat["GWGC"], cat["PGC"], cat["HyperLEDA"]
        _2MASS, SDSS = cat["_2MASS"], cat["SDSS-DR12"]

        with h5py.File(catalogFile, 'w') as f:
            f.create_dataset('ra', data=ra)
            f.create_dataset('dec', data=dec)
            f.create_dataset('distmpc', data=distmpc)
            f.create_dataset('Flag1', data=flag1)
            f.create_dataset('magb', data=magb)
            f.create_dataset('BMAG', data=BMAG)
            f.create_dataset('Jmag', data=Jmag)
            f.create_dataset('Hmag', data=Hmag)
            f.create_dataset('Kmag', data=Kmag)
            f.create_dataset('Flag2', data=flag2)
            f.create_dataset('Flag3', data=flag3)

            f.create_dataset('z', data=z)
            # Add galaxy identifier
            f.create_dataset('GWGC', data=GWGC)
            f.create_dataset('PGC', data=PGC)
            f.create_dataset('HyperLEDA', data=HyperLEDA)
            f.create_dataset('2MASS', data=_2MASS)
            f.create_dataset('SDSS', data=SDSS)


def glade_query(ra_deg, dec_deg, rad_deg, dist_constraint=[0,3000], maxmag=20,
               maxsources=-1, online=True, catalogFile='catalogs/glade.hdf5'):
    """
    Query glade @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                dist_constraint: min and max distance in Mpc
                maxmag: upper limit magnitude (optional)
                maxsources: maximum number of sources
                online: True: online query - False: use catalog locally
                catalogFile: path to glade catalog

    Use online because offline conesearch returns different results
    returns: astropy.table object
    """
    if online:

        vquery = Vizier(columns=['PGC', 'GWGC', 'HyperLEDA', '2MASS', 'SDSS-DR12',
                                 'Flag1', 'RAJ2000', 'DEJ2000', 'Dist', 'z',
                                 'Bmag', 'BMAG', 'Jmag', 'Hmag', 'Kmag',
                                 'Flag2' , 'Flag3'],
                        column_filters={"Dist": ">%f & <%f" % (dist_constraint[0], dist_constraint[1])},
                        row_limit = maxsources)
        field = SkyCoord(ra=ra_deg, dec=dec_deg,
                         unit=(u.deg, u.deg),
                         frame='icrs')
        table = vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="VII/281/glade2")[0]
    else:
        from astropy.table import Table

        with h5py.File(catalogFile, 'r') as f:
            ra, dec = f['ra'][:], f['dec'][:]
            distmpc, z = f['distmpc'][:], f['z'][:]
            magb,  BMAG = f['magb'][:],  f['BMAG'][:]
            Jmag, Hmag= f['Jmag'][:], f['Hmag'][:]
            Kmag = f['Kmag'][:]
            Flag1, Flag2, Flag3 = f['Flag1'][:], f['Flag2'][:], f['Flag3'][:]
            GWGC, PGC, _2MASS = f['GWGC'][:], f['PGC'][:], f['2MASS'][:]
            HyperLEDA, SDSS = f['HyperLEDA'][:], f['SDSS'][:]
            # Convert bytestring to unicode
            GWGC = GWGC.astype('U')
            PGC = PGC.astype('U')
            HyperLEDA = HyperLEDA.astype('U')
            _2MASS = _2MASS.astype('U')
            SDSS = SDSS.astype('U')

        table = Table([PGC, GWGC, HyperLEDA, _2MASS, SDSS, Flag1, ra, dec, distmpc,
                       z, magb, BMAG, Jmag, Hmag, Kmag, Flag2, Flag3],
                      names = ['PGC', 'GWGC','HyperLEDA','_2MASS', 'SDSS-DR12','Flag1','RAJ2000',
                              'DEJ2000', 'Dist', 'z', 'Bmag', 'BMAG',
                              'Jmag', 'Hmag', 'Kmag', 'Flag2', 'Flag3'])
        # constraint on distance
        mask = (table['Dist'] > dist_constraint[0]) & (table['Dist'] < dist_constraint[1])
        table = table[mask]
        # cone search
        mask = (table['RAJ2000'] - ra_deg)**2 + (table['DEJ2000'] - dec_deg)**2 <= rad_deg**2
        table=table[mask]

    return table
