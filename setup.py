from setuptools import setup
import rnaseqtools

setup(name='rnaseqtools',
      version=0.1,
      description='Handling RNAseq data, wrappers for aligners, functions on bamfiles',
      url='http://github.com/redst4r/rnaseqtools/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='RNAseq',
      packages=['rnaseqtools'],
      install_requires=[
          'toolz',
          'numpy',
          'pysam',
          'tqdm',
          'scipy',
          'pandas',
          'bioservices',
          'pyliftover',
          ],
      zip_safe=False)
