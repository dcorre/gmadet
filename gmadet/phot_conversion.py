#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def poly(x, coefficients):
    """
     Compute a polynome, useful for transformation laws

     parameters:
         x: Variable for the polynome
         coefficients: list of coefficients

     returns: float result of the polynome
     """
    poly = 0
    for i, coef in enumerate(coefficients):
        poly += coef * x ** i

    return poly


def gaia2Johnson(band, Gaia_Table):
    """
    Give the transformation laws to go from gaia photometric system
    to Johnson-Cousins photometric system
    Transformation given by https://gea.esac.esa.int/archive/documentation/
    GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/
    ssec_cu5pho_PhotTransf.html
    Table 5.8

    parameters: band: filter used to get the image
                Gaia_Table: astropy.table with information from Gaia for
                the sources in the image

    returns: astropy.table with informations from Gaia and the computation
             for the filter of the image
    """
    V_band = ["V"]  # Differen15.690t names that can be
    R_band = ["R"]  # found in header
    I_band = ["I"]
    B_band = ["B"]

    if band in V_band:
        # Validity domain
        mask = (-0.5 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 2.75)
        Result = Gaia_Table[mask]
        coefficients = [-0.01760, -0.006860, -0.1732]
        Result["VMag"] = Result["phot_g_mean_mag"] - \
            poly(Result["bp_rp"], coefficients)
        Result["calib_err"] = 0.045858

    if band in R_band:
        # Validity domain
        mask = (-0.5 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 2.75)
        Result = Gaia_Table[mask]
        coefficients = [-0.003226, 0.3833, -0.1345]
        Result["RMag"] = Result["phot_g_mean_mag"] - \
            poly(Result["bp_rp"], coefficients)
        Result["calib_err"] = 0.04840

    if band in I_band:
        # Validity domain
        mask = (-0.5 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 2.75)
        Result = Gaia_Table[mask]
        coefficients = [0.02085, 0.7419, -0.09631]
        Result["IMag"] = Result["phot_g_mean_mag"] - poly(
            Result["bp_rp"], coefficients
        )  # Revoir cette trabsformation, pt etre mal  implementee
        Result["calib_err"] = 0.04956

    if band in B_band:
        raise ValueError("No transformation law for Gaia->B_Jonhson")

    return Result


def gaia2SDSS(band, Gaia_Table):
    """
    Give the transformation laws to go from gaia photometric system
    to SDSS photometric system
    Transformation given by https://gea.esac.esa.int/archive/documentation/
    GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/
    ssec_cu5pho_PhotTransf.html
    Table 5.7

    parameters: band: filter used to get the image
                Gaia_Table: astropy.table with information from Gaia for
                the sources in the image

    returns: astropy.table with informations from Gaia and the computation
             for the filter of the image
    """
    g_band = ["g", "g'", "sdss g", "SDSS g"]  # Different names that can be
    r_band = ["r", "r'", "sdss r", "SDSS r"]  # found in header
    i_band = ["i", "i'", "sdss i", "SDSS i"]
    z_band = ["z", "z'", "sdss z", "SDSS z"]

    if band in r_band:
        # Validity domain
        mask = (0.2 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 2.7)
        Result = Gaia_Table[mask]
        coefficients = [-0.12879, 0.24662, -0.027464, -0.049465]
        Result["r_SDSSMag"] = Result["phot_g_mean_mag"] - poly(
            Result["bp_rp"], coefficients
        )
        Result["calib_err"] = 0.066739

    if band in i_band:
        # Validity domain
        mask = (0.0 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 4.5)
        Result = Gaia_Table[mask]
        coefficients = [-0.29676, 0.64728, -0.10141, 0.0]
        Result["i_SDSSMag"] = Result["phot_g_mean_mag"] - poly(
            Result["bp_rp"], coefficients
        )
        Result["calib_err"] = 0.098957

    if band in g_band:
        # Validity domain
        mask = (-0.5 < Gaia_Table["bp_rp"]) & (Gaia_Table["bp_rp"] < 2.0)
        Result = Gaia_Table[mask]
        coefficients = [0.13518, -0.46245, -0.25171, 0.021349]
        Result["g_SDSSMag"] = Result["phot_g_mean_mag"] - poly(
            Result["bp_rp"], coefficients
        )
        Result["calib_err"] = 0.16497

    if band in z_band:
        raise ValueError("No transformation law for Gaia->z_sdss")

    return Result


def usno2Johnson(band, USNO_Table):
    """
    Give the transformation laws to go from USNO photometric system
    to Johnson photometric system
    Transformation given by http://www.mpe.mpg.de/~jcg/GROND/calibration.html

    parameters: band: filter used to get the image
                USNO_Table: astropy.table with information from USNO for
                the sources in the image

    returns: astropy.table with informations from USNO and the computation
             for the filter of the image
    """
    V_band = ["V"]  # Different names that can be
    R_band = ["R"]  # found in header
    I_band = ["I"]
    B_band = ["B"]

    if band in V_band:
        USNO_Table["V_JohnsonMag"] = (
            0.444 * USNO_Table["B1mag"] + 0.556 * USNO_Table["R1mag"]
        )
        USNO_Table["calib_err"] = 0.5
    if band in B_band:
        USNO_Table["B_JohnsonMag"] = USNO_Table["B1mag"]
        USNO_Table["calib_err"] = 0.5
    if band in R_band:
        USNO_Table["R_JohnsonMag"] = USNO_Table["R1mag"]
        USNO_Table["calib_err"] = 0.5
    if band in I_band:
        USNO_Table["I_JohnsonMag"] = USNO_Table["Imag"]
        USNO_Table["calib_err"] = 0.5

    return USNO_Table


def SDSS2Johnson(band, SDSS_Table):
    """
    Give the transformation laws to go from SDSS photometric system
    to Johnson photometric system
    Transformation given by http://www.sdss3.org/dr8/algorithms/
    sdssUBVRITransform.php

    parameters: band: filter used to get the image
                SDSS_Table: astropy.table with information from SDSS for
                the sources in the image

    returns: astropy.table with informations from SDSS and the computation
             for the filter of the image
    """
    V_band = ["V"]  # Different names that can be
    R_band = ["R"]  # found in header
    I_band = ["I"]
    B_band = ["B"]

    if band in V_band:
        coefficients = [-0.016, -0.573]
        SDSS_Table["g-r"] = SDSS_Table["gmag"] - SDSS_Table["rmag"]
        SDSS_Table["VMag"] = SDSS_Table["gmag"] + \
            poly(SDSS_Table["g-r"], coefficients)
        SDSS_Table["calib_err"] = (
            np.sqrt(
                (0.002 / 0.573) ** 2
                + (0.002 / 0.016) ** 2
                + (SDSS_Table["e_gmag"] / SDSS_Table["gmag"]) ** 2
                + (SDSS_Table["e_rmag"] / SDSS_Table["rmag"]) ** 2
            )
            * SDSS_Table["VMag"]
        )

    if band in R_band:
        coefficients = [0.152, -0.257]
        SDSS_Table["r-i"] = SDSS_Table["rmag"] - SDSS_Table["imag"]
        SDSS_Table["RMag"] = SDSS_Table["rmag"] + \
            poly(SDSS_Table["r-i"], coefficients)
        SDSS_Table["calib_err"] = (
            np.sqrt(
                (0.004 / 0.257) ** 2
                + (0.002 / 0.152) ** 2
                + (SDSS_Table["e_rmag"] / SDSS_Table["rmag"]) ** 2
                + (SDSS_Table["e_imag"] / SDSS_Table["imag"]) ** 2
            )
            * SDSS_Table["RMag"]
        )

    if band in I_band:
        coefficients = [-0.394, -0.409]
        SDSS_Table["i-z"] = SDSS_Table["imag"] - SDSS_Table["zmag"]
        SDSS_Table["IMag"] = SDSS_Table["imag"] + \
            poly(SDSS_Table["i-z"], coefficients)
        SDSS_Table["calib_err"] = (
            np.sqrt(
                (0.006 / 0.409) ** 2
                + (0.002 / 0.394) ** 2
                + (SDSS_Table["e_zmag"] / SDSS_Table["zmag"]) ** 2
                + (SDSS_Table["e_imag"] / SDSS_Table["imag"]) ** 2
            )
            * SDSS_Table["IMag"]
        )

    if band in B_band:
        coefficients = [0.219, 0.312]
        SDSS_Table["g-r"] = SDSS_Table["gmag"] - SDSS_Table["rmag"]
        SDSS_Table["BMag"] = SDSS_Table["gmag"] + \
            poly(SDSS_Table["g-r"], coefficients)
        SDSS_Table["calib_err"] = (
            np.sqrt(
                (0.002 / 0.219) ** 2
                + (0.003 / 0.312) ** 2
                + (SDSS_Table["e_gmag"] / SDSS_Table["gmag"]) ** 2
                + (SDSS_Table["e_rmag"] / SDSS_Table["rmag"]) ** 2
            )
            * SDSS_Table["BMag"]
        )

    return SDSS_Table


def PS2Johnson(band, PS_Table):
    """
    Give the transformation laws to go from Pan-STARRS photometric system
    to Johnson photometric system
    Transformation given by http://www.sdss3.org/dr8/algorithms/
    sdssUBVRITransform.php

    parameters: band: filter used to get the image
                PS_Table: astropy.table with information from PS for
                the sources in the image

    returns: astropy.table with informations from PS and the computation
             for the filter of the image
    """
    V_band = ["V"]  # Different names that can be
    R_band = ["R"]  # found in header
    I_band = ["I"]
    B_band = ["B"]

    if band in V_band:
        coefficients = [-0.016, -0.573]
        PS_Table["g-r"] = PS_Table["gmag"] - PS_Table["rmag"]
        PS_Table["VMag"] = PS_Table["gmag"] + \
            poly(PS_Table["g-r"], coefficients)
        PS_Table["calib_err"] = (
            np.sqrt(
                (0.002 / 0.573) ** 2
                + (0.002 / 0.016) ** 2
                + (PS_Table["e_gmag"] / PS_Table["gmag"]) ** 2
                + (PS_Table["e_rmag"] / PS_Table["rmag"]) ** 2
            )
            * PS_Table["VMag"]
        )

    if band in R_band:
        coefficients = [0.152, -0.257]
        PS_Table["r-i"] = PS_Table["rmag"] - PS_Table["imag"]
        PS_Table["RMag"] = PS_Table["rmag"] + \
            poly(PS_Table["r-i"], coefficients)
        PS_Table["calib_err"] = (
            np.sqrt(
                (0.004 / 0.257) ** 2
                + (0.002 / 0.152) ** 2
                + (PS_Table["e_rmag"] / PS_Table["rmag"]) ** 2
                + (PS_Table["e_imag"] / PS_Table["imag"]) ** 2
            )
            * PS_Table["RMag"]
        )

    if band in I_band:
        coefficients = [-0.394, -0.409]
        PS_Table["i-z"] = PS_Table["imag"] - PS_Table["zmag"]
        PS_Table["IMag"] = PS_Table["imag"] + \
            poly(PS_Table["i-z"], coefficients)
        PS_Table["calib_err"] = (
            np.sqrt(
                (0.006 / 0.409) ** 2
                + (0.002 / 0.394) ** 2
                + (PS_Table["e_zmag"] / PS_Table["zmag"]) ** 2
                + (PS_Table["e_imag"] / PS_Table["imag"]) ** 2
            )
            * PS_Table["IMag"]
        )

    if band in B_band:
        coefficients = [0.219, 0.312]
        PS_Table["g-r"] = PS_Table["gmag"] - PS_Table["rmag"]
        PS_Table["BMag"] = PS_Table["gmag"] + \
            poly(PS_Table["g-r"], coefficients)
        PS_Table["calib_err"] = (
            np.sqrt(
                (0.002 / 0.219) ** 2
                + (0.003 / 0.312) ** 2
                + (PS_Table["e_gmag"] / PS_Table["gmag"]) ** 2
                + (PS_Table["e_rmag"] / PS_Table["rmag"]) ** 2
            )
            * PS_Table["BMag"]
        )

    return PS_Table
