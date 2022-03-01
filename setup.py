#!/usr/bin/env python3


"""Py-SP2

A package for analyzing Single Particle Soot Photometer data

"""


DOCLINES = __doc__.split("\n")

import glob

from setuptools import setup, find_packages
from os import path

# Needed to build Windows distribution on Linux
# Work around mbcs bug in distutils.
# http://bugs.python.org/issue10945
import codecs
try:
    codecs.lookup('mbcs')
except LookupError:
    ascii = codecs.lookup('ascii')
    func = lambda name, enc=ascii: {True: enc}.get(name=='mbcs')
    codecs.register(func)

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

NAME = 'pysp2'
MAINTAINER = 'Bobby Jackson'
MAINTAINER_EMAIL = 'rjackson@anl.gov'
URL = ''
DESCRIPTION = DOCLINES[0]
LONG_DESCRIPTION = "\n".join(DOCLINES[2:])
LICENSE = 'BSD'
PLATFORMS = "Linux, Windows, OSX"
MAJOR = 1
MINOR = 0
MICRO = 0
#SCRIPTS = glob.glob('scripts/*')
#TEST_SUITE = 'nose.collector'
#TESTS_REQUIRE = ['nose']
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

def setup_package():
    """ Setup of PySP2  package. """
    setup(
        name=NAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        url=URL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type='text/markdown',
        version=VERSION,
        license=LICENSE,
        platforms=PLATFORMS,
        include_package_data=True,
        install_requires=['matplotlib', 'pandas', 'numpy', 'scipy'],
        packages=find_packages(exclude=['contrib', 'docs', 
                                        'tests', 'examples']),
        project_urls={
            'Source': 'https://github.com/rcjackson/PySP2'},
        #test_suite=TEST_SUITE,
        #tests_require=TESTS_REQUIRE,
        #scripts=SCRIPTS,
    )

if __name__ == '__main__':
    setup_package()
