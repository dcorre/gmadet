#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gmadet` package."""

import pytest
from gmadet.cli.psf import main as gmadet_psf


def test_psf(script_runner):
    """test gmadet-psf"""
    ret = script_runner.run('gmadet-psf', '--path_data',
            'gmadet/data_test/ATLAS18qqn-S001-R001-C001-SDSS_g.fits',
            '--telescope', 'IRIS')
    assert ret.success
    assert ret.stderr == ''
    
