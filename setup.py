#!/usr/bin/env python
# coding=utf-8

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

setup(
    name='hydrotopo',
    version='0.1.2',
    description='A project which find topologic relations among stations along rivers',
    entry_points={'console_scripts': ['calcstream = hydrotopo.cli:main']},
    long_description='https://readthedocs.org/projects/calctest/',
    author='forestbat',
    url='https://github.com/iHeadWater/HydroTopo',
    packages=['hydrotopo'],
    package_dir={'hydrotopo': 'hydrotopo'},
    include_package_data=True,
    install_requires=['geopandas>=0.12.2', 'pytest>=7.1.3', 'setuptools>=70.0.0', 'click>=8.1.3', 'pandas>=1.4.3',
                      'shapely>=2.0.1', 'pyogrio>=0.4.2', 'igraph>=0.10.4'],
    license='MIT',
    zip_safe=False,
    keywords='hydrotopo',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)