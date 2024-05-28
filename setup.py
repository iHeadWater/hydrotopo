#!/usr/bin/env python
# coding=utf-8

import os
import sys
import random
import string

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

setup(
    name='hydro_topo',
    version='0.0.2',
    description='A project which find topologic relations among stations along rivers',
    entry_points={'console_scripts': ['calcstream = hydro_topo.cli:main']},
    long_description='https://readthedocs.org/projects/calctest/',
    author='forestbat',
    url='https://github.com/iHeadWater/HydroTopo',
    packages=['hydro_topo'],
    package_dir={'hydro_topo': 'hydro_topo'},
    include_package_data=True,
    install_requires=['geopandas>=0.12.2', 'pytest>=7.1.3', 'setuptools>=63.4.1', 'click>=8.1.3', 'pandas>=1.4.3',
                      'shapely>=2.0.1', 'pyogrio>=0.4.2', 'igraph>=0.10.4'],
    license='MIT',
    zip_safe=False,
    keywords='hydro_topo',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
