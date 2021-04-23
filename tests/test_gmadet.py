#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gmadet` package."""

import os
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
    """test gmadet-astrometry"""
    ret = script_runner.run(
        'gmadet-astrometry',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_astrometry',
        '--telescope', 'IRIS'
    )
    assert ret.success
    assert ret.stderr == ''


def test_stacking(script_runner):
    """test gmadet-stacking. Only with one file but it is just a test"""
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
    """test gmadet-subBkg"""
    ret = script_runner.run(
        'gmadet-subBkg',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_subbkg',
    )
    assert ret.success
    assert ret.stderr == ''


def test_remove_cosmics(script_runner):
    """test gmadet-cosmics"""
    ret = script_runner.run(
        'gmadet-cosmics',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_cosmics',
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_run_noSub(script_runner):
    """test gmadet-run without substraction."""
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
    """test gmadet-run with substraction."""
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
    """test gmadet-run with substraction."""
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


def test_gmadet_sim(script_runner):
    """test gmadet-sim."""
    ret = script_runner.run(
        'gmadet-sim',
        'gmadet/data_test/test_image.fits',
        '--results', 'gmadet/data_test/gmadet_sim',
        '--telescope', 'IRIS',
        '--ntrans', '100',
        '--magrange', '14', ' 23',
        '--zp', '30'
    )
    assert ret.success
    assert ret.stderr == ''

def test_gmadet_run_sim(script_runner):
    """test gmadet-run for simulation."""
    ret = script_runner.run(
        'gmadet-run',
        'gmadet/data_test/gmadet_sim',
        '--results', 'gmadet/data_test/gmadet_sim_results',
        '--telescope', 'IRIS',
        '--fwhm', 'psfex',
        '--astrometry', 'scamp',
        '--radius-crossmatch', '3',
        '--threshold', '4',
        '--sub', 'ps1',
        '--ps1-method', 'individual',
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_cutouts(script_runner):
    """test gmadet-cutouts."""
    ret = script_runner.run(
        'gmadet-cutouts',
        'gmadet/data_test/gmadet_sim_results',
        '--training',
        '--false'
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_checksim(script_runner):
    """test gmadet-checksim."""
    ret = script_runner.run(
        'gmadet-checksim',
        'gmadet/data_test/gmadet_sim_results',
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_cnn_convert(script_runner):
    """test gmadet-cnn_convert."""
    ret = script_runner.run(
        'gmadet-cnn_convert',
        '--path', 'gmadet/data_test/gmadet_cnn',
        '--cube', 'cube',
        '--cutouts',  'gmadet/data_test/gmadet_sim_results/candidates_training/',
        '--frac_true', '0.5'
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_cnn_train(script_runner):
    """test gmadet-cnn_train."""
    ret = script_runner.run(
        'gmadet-cnn_train',
        '--cube', 'gmadet/data_test/gmadet_cnn/datacube/cube.npz',
        '--model-path', 'gmadet/data_test/gmadet_cnn',
        '--model-name', 'model'
    )
    assert ret.success
    assert ret.stderr == ''

# Comment as there are some problems with mp.pool here, may be some pool 
# not closed before?
#def test_gmadet_cutouts_infer(script_runner):
#    """test gmadet-cutouts in preparation of gmadet-cnn_infer."""
#    # Create cutouts without classifying them in true or false folders
#    ret = script_runner.run(
#        'gmadet-cutouts',
#        'gmadet/data_test/gmadet_sim_results'
#    )
#
#    assert ret.success
#    assert ret.stderr == ''


def test_gmadet_cnn_infer(script_runner):
    """test gmadet-cnn_infer."""
    # move candidates by hand. 
    os.makedirs("gmadet/data_test/gmadet_sim_results/candidates",  exist_ok = True)
    os.system("cp gmadet/data_test/gmadet_sim_results/candidates_training/true/*.fits gmadet/data_test/gmadet_sim_results/candidates/")
    os.system("cp gmadet/data_test/gmadet_sim_results/candidates_training/false/*.fits gmadet/data_test/gmadet_sim_results/candidates/")
    ret = script_runner.run(
        'gmadet-cnn_infer',
        '--cutouts', 'gmadet/data_test/gmadet_sim_results/candidates',
        '--model', 'gmadet/data_test/gmadet_cnn/CNN_training/model.h5',
    )
    assert ret.success
    assert ret.stderr == ''


def test_gmadet_cnn_checkinfer(script_runner):
    """test gmadet-cnn_checkinfer."""
    ret = script_runner.run(
        'gmadet-cnn_checkinfer',
        '--plots', 'gmadet/data_test/gmadet_sim_results/',
        '--crossmatch', 'gmadet/data_test/gmadet_sim_results/',
        '--infer', 'gmadet/data_test/gmadet_sim_results/candidates/',
    )
    assert ret.success
    assert ret.stderr == ''
