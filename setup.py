#!/usr/bin/env python2.7

from distutils.core import setup

setup(
    name='snppipeline',
    version='0.1.1',
    description='Script and functions for SNP matrix construction',
    author='Yan Luo',
    author_email='yan.luo@fda.hhs.gov',
    #url='http://', #TODO
    packages=['snppipeline'],
    long_description="""
    'snppipeline' is a pipeline for the production of SNP matrices for the
    phylogenetic analysis of pathogenic organisms sequenced from various
    samples of interest to food safety.
    """,
    #download_url='',#TODE    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords='networking eventlet nonblocking internet',
    license='GPL',
    install_requires=[
        'setuptools',
    ],
)
