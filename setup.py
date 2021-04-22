#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
        # tensorflow 2.3.2 requires numpy<1.19.0,>=1.16.0
        'numpy==1.18.5',
        'astropy',
        'matplotlib',
        'pandas',
        'astroquery',
        'astroML',
        'cython',
        'Sphinx',
        'twine',
        'scipy',
        'pandas',
        'shapely',
        'requests',
        'h5py',
        'scikit-image',
        'lacosmic',
        'hjson',
        'voevent-parse',
        'xmltodict',
        'astroML',
        'photutils',
        'keras',
        'keras-vis',
        'tensorflow==2.3.2',
        'regions',
        'opencv-python-headless',
        'astroquery',
        'astroscrappy'
       ]
setup_requirements = [
        'pytest-runner',
        'flake8',
        'bumpversion',
        'wheel',
        'twine']

test_requirements = [
        'pytest',
        'pytest-cov',
        'pytest-console-scripts',
        'pytest-html',
        'watchdog']

setup(
    author="David Corre",
    author_email='david.corre.fr@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Tools to help identification of transients "
                "for Grandma network",
    entry_points={
        'console_scripts': [
            'gmadet-run = gmadet.cli.run_gmadet:main',
            'gmadet-astrometry = gmadet.cli.astrometry:main',
            'gmadet-stacking = gmadet.cli.stacking:main',
            'gmadet-psf = gmadet.cli.psf:main',
            'gmadet-subBkg = gmadet.cli.subBkg:main',
            'gmadet-cosmics = gmadet.cli.cosmics:main',
            'gmadet-sim = gmadet.cli.sim:main',
            'gmadet-cutouts = gmadet.cli.make_cutouts:main',
            'gmadet-checksim = gmadet.cli.checksim:main',
            'gmadet-cnn_convert = gmadet.cli.cnn_convert:main',
            'gmadet-cnn_train = gmadet.cli.cnn_train:main',
            'gmadet-cnn_infer = gmadet.cli.cnn_infer:main',
            'gmadet-cnn_checkinfer = gmadet.cli.cnn_checkinfer:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords=['gmadet', 'transients', 'detection pipeline', 'astronomy',
              'image substraction', 'CNN'],
    name='gmadet',
    # packages=find_packages(include=['gmadet']),
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dcorre/gmadet',
    version='0.1.0',
    zip_safe=False,
)
