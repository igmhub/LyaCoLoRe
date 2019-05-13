#!/usr/bin/env python

import glob

from setuptools import setup

scripts = glob.glob('bin/*')

description = "Fast generation of Lya catalogs"

version="0.1"
setup(name="LyaCoLoRe",
    version=version,
    description=description,
    url="https://github.com/igmhub/LyaCoLoRe",
    author="<fill the names there>",
    author_email="<fill the email there>",
    packages=['LyaCoLoRe'],
    package_dir = {'': 'py'},
    install_requires=['numpy','scipy','iminuit','healpy','fitsio',
        'numba','future','setuptools'],
    test_suite='LyaCoLoRe.test.test_cor',
    scripts = scripts
)
