#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""Main module."""

from sources_det import sources_det
from crossmatch import crossmatch

# Define path to data
path_data = '/home/corre/codes/algoCNN/rawdata/OAJ-T80/'

# Run sextractor to find sources in all images
sources_det(path_data)

# Cross match with catalogs
#crossmatch('sources_cat/', catalog='PS1')
