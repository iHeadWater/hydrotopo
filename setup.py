#!/usr/bin/env python
#coding=utf-8

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
    name='dijkstra-conda',
    version=''.join(random.sample(string.digits, 8)),
    description='A test project',
    entry_points={'console_scripts': ['calcstream = dijkstra-conda.cli:main']},
    long_description='https://readthedocs.org/projects/calctest/',
    author='dijkstra-conda-test',
    url='https://github.com/iHeadWater/dijkstra-conda',
    packages=[
        'dijkstra-conda',
    ],
    package_dir={'dijkstra-conda': 'dijkstra-conda'},
    include_package_data=True,
    install_requires=[
    ],
    license='MIT',
    zip_safe=False,
    keywords='dijkstra-conda',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
