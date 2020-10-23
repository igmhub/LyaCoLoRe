from distutils.core import setup
import glob

scripts = glob.glob('scripts/*')

setup(name='LyaCoLoRe',
      version='0.1',
      description='Fast generation of Lya catalogs',
      author='igmhub',
      author_email='github.com/igmhub',
      url='https://github.com/igmhub/LyaCoLoRe',
      packages=['lyacolore'],
      package_dir = {'lyacolore': 'py'},
      test_suite='picca.test.test_cor',
      scripts = scripts
     )
