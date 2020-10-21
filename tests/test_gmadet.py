#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gmadet` package."""

import pytest
from gmadet.cli.psf import main as gmadet_psf


def test_psf(script_runner):
    """test gmadet-psf"""
    ret = script_runner.run(
        'gmadet-psf',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_psf',
        '--telescope', 'IRIS'
    )
    assert ret.success
    assert ret.stderr == ''


def test_astrometry(script_runner):
    "test gmadet-astrometry"
    ret = script_runner.run(
        'gmadet-astrometry',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_astrometry',
        '--telescope', 'IRIS'
    )
    assert ret.success
    assert ret.stderr == ''


def test_stacking(script_runner):
    "test gmadet-stacking. Only with one file but it is just a test"
    ret = script_runner.run(
        'gmadet-stacking',
        'gmadet/data_test/stacking/',
        '--results', 'gmadet/data_test/gmadet_stacking',
        '--radius', '1',
        '--deltat', '0.5'
    )
    assert ret.success
    assert ret.stderr == ''


def test_subBkg(script_runner):
    "test gmadet-subBkg"
    ret = script_runner.run(
        'gmadet-subBkg',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_subbkg',
    )
    assert ret.success
    assert ret.stderr == ''


def test_remove_cosmics(script_runner):
    "test gmadet-cosmics"
    ret = script_runner.run(
        'gmadet-cosmics',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_cosmics',
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_run_noSub(script_runner):
    "test gmadet-run without substraction."
    ret = script_runner.run(
        'gmadet-run',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_results',
        '--telescope', 'IRIS',
        '--fwhm', 'psfex',
        '--astrometry', 'scamp',
        '--radius-crossmatch', '3',
        '--threshold', '4',
    )
    assert ret.success
    # Tensorflow is displaying a lot of warning that make the 
    # line below crashing. Do not understand why only this
    # test though...
    # assert ret.stderr == ''


def test_gmadet_run_Sub_individual(script_runner):
    "test gmadet-run with substraction."
    ret = script_runner.run(
        'gmadet-run',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_results',
        '--telescope', 'IRIS',
        '--fwhm', 'psfex',
        '--astrometry', 'scamp',
        '--radius-crossmatch', '3',
        '--threshold', '4',
        '--sub', 'ps1',
        '--ps1-method', 'individual',
        '--mosaic'
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_run_Sub_mosaic(script_runner):
    "test gmadet-run with substraction."
    ret = script_runner.run(
        'gmadet-run',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_results',
        '--telescope', 'IRIS',
        '--fwhm', 'psfex',
        '--astrometry', 'scamp',
        '--radius-crossmatch', '3',
        '--threshold', '4',
        '--sub', 'ps1',
        '--ps1-method', 'mosaic',
    )
    assert ret.success
    assert ret.stderr == ''
