#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gmadet` package."""

import pytest
from gmadet.cli.psf import main as gmadet_psf


def test_psf(script_runner):
    """test gmadet-psf"""
    ret = script_runner.run(
            'gmadet-psf',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            '--telescope', 'IRIS'
            )
    assert ret.success
    assert ret.stderr == ''

def test_astrometry(script_runner):
    "test gmadet-astrometry"
    ret = script_runner.run(
            'gmadet-astrometry',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            '--telescope', 'IRIS'
            )
    assert ret.success
    assert ret.stderr == ''

def test_stacking(script_runner):
    "test gmadet-stacking. Only with one file but it is just a test"
    ret = script_runner.run(
            'gmadet-stacking',
            '--path_data', 'gmadet/data_test/stacking/',
            '--radius', '1',
            '--deltaT', '0.5'
            )
    assert ret.success
    assert ret.stderr == ''

def test_subBkg(script_runner):
    "test gmadet-subBkg"
    ret = script_runner.run(
            'gmadet-subBkg',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            )
    assert ret.success
    assert ret.stderr == ''

def test_remove_cosmics(script_runner):
    "test gmadet-cosmics"
    ret = script_runner.run(
            'gmadet-cosmics',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            )
    assert ret.success
    assert ret.stderr == ''

def test_gmadet_run_noSub(script_runner):
    "test gmadet-run without substraction."
    ret = script_runner.run(
            'gmadet-run',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            '--telescope', 'IRIS',
            '--FWHM', 'psfex',
            '--doAstrometry', 'scamp',
            '--radius_crossmatch', '3',
            '--threshold', '4',
            )
    assert ret.success
    assert ret.stderr == ''

def test_gmadet_run_Sub_individual(script_runner):
    "test gmadet-run with substraction."
    ret = script_runner.run(
            'gmadet-run',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            '--telescope', 'IRIS',
            '--FWHM', 'psfex',
            '--doAstrometry', 'scamp',
            '--radius_crossmatch', '3',
            '--threshold', '4',
            '--doSub', 'ps1',
            '--ps1_method', 'individual', 
            '--doMosaic'
            )
    assert ret.success
    assert ret.stderr == ''

def test_gmadet_run_Sub_mosaic(script_runner):
    "test gmadet-run with substraction."
    ret = script_runner.run(
            'gmadet-run',
            '--path_data',
            'gmadet/data_test/test_image.fits',
            '--telescope', 'IRIS',
            '--FWHM', 'psfex',
            '--doAstrometry', 'scamp',
            '--radius_crossmatch', '3',
            '--threshold', '4',
            '--doSub', 'ps1',
            '--ps1_method', 'mosaic',
            )
    assert ret.success
    assert ret.stderr == ''
