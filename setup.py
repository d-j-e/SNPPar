#!/usr/bin/env python

from distutils.core import setup

setup(
    name='snppar',
    version='0.4.2dev',
    author='David Edwards',
    author_email='David.Edwards@monash.edu',
    packages=['snppar'],
    scripts=['scripts/snppar.py'],
    entry_points={
        'console_scripts': ['snppar = snppar:main']
    },
    package_dir = {'snppar': 'scripts'},
    url='https://github.com/d-j-e/SNPPar',
    license='LICENSE.txt',
    description='Parallel/Homoplasic SNP Finder for Bacteria',
    long_description=('This program is designed to take a SNP'
                      'alignment and a phylogenetic tree in order to'
                      'find any homoplasic SNPs, and define them if required)'
                      'By default, reports both the homoplasic events and '
                      'all mutation events. Also maps the SNPs to the tree.'),
    install_requires=[
        'biopython>=1.66',
        'ete3',
        'phylo-treetime'
    ],
)
