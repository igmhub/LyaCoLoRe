#!/usr/bin/env python

from setuptools import setup, find_packages
import glob

scripts = glob.glob('bin/*')

requirements = ['numpy','scipy','iminuit','healpy','fitsio',
    'numba','future','setuptools','configargparse']

setup(name='LyaCoLoRe',
      version='0.1',
      description='Fast generation of Lya catalogs',
      author='igmhub',
      author_email='james.farr.17@ucl.ac.uk',
      url='https://github.com/igmhub/LyaCoLoRe',
      packages=find_packages('py'),
      package_dir = {'':'py'},
      install_requires=requirements,
      test_suite='lyacolore.test.test_cor',
      scripts = scripts
     )
