#!/usr/bin/env python

import glob

from setuptools import setup

scripts = glob.glob('bin/*')

description = "Fast generation of Lya catalogs"

version='0.1'
setup(name="LyaCoLoRe",
    version=version,
    description=description,
    url="https://github.com/igmhub/LyaCoLoRe",
    author="<your name here>",
    author_email="<your email here>",
    packages=['LyaCoLoRe'],
    package_dir = {'': 'py'},
    package_data = {'LyaCoLoRe': ['etc/']},
    install_requires=['numpy','scipy','iminuit','healpy','fitsio',
        'numba','future','setuptools'],
    test_suite='LyaCoLoRe.test.test_cor',
    scripts = scripts
    )
