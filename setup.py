#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import with_statement

import io
import sys
import warnings
from setuptools import setup
from setuptools import find_packages

VERSION = 'undefined'
install_requires = [
    'six',
    'intervaltree>=2.1.0',
    'pyfasta>=0.5.2',
    'numpy>=1.11.1',
]

extra = {}

with open('veppy/version.py') as f:
    for row in f.readlines():
        if row.startswith('VERSION'):
            exec(row)

if sys.version_info < (2, 6):
    warnings.warn(
        'Python 2.5 is not officially supported by veppy. '
        'If you have any questions, please file an issue on GitHub.',
        DeprecationWarning)
elif sys.version_info >= (3, 0):
    extra['use_2to3'] = True

with io.open('README.md', encoding="utf-8") as f:
    long_description = f.read()

setup(
    name='veppy',
    version=VERSION,
    description='Genetic Variant Effect Prediction for Python',
    long_description=long_description,
    author='Solve, Inc.',
    author_email='contact@solvebio.com',
    url='https://github.com/solvebio/veppy',
    include_package_data=True,
    packages=find_packages(),
    package_dir={"veppy": "veppy"},
    platforms='any',
    install_requires=install_requires,
    # TODO: Test suite
    # test_suite='veppy.test.all',
    entry_points={
        # TODO: Command-line veppy
        # 'console_scripts': ['veppy = veppy.cli.main:main']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    **extra
)
