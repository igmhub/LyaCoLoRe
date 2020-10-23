from setuptools import setup
import glob

scripts = glob.glob('scripts/*')

setup(name='LyaCoLoRe',
      version='0.1',
      description='Fast generation of Lya catalogs',
      author='igmhub',
      author_email='james.farr.17@ucl.ac.uk',
      url='https://github.com/igmhub/LyaCoLoRe',
      packages=['lyacolore'],
      package_dir = {'lyacolore': 'py'},
      install_requires=['numpy','scipy','iminuit','healpy','fitsio',
          'numba','future','setuptools','configargparse'],
      test_suite='picca.test.test_cor',
      scripts = scripts
     )
