#!/usr/bin/env python2.7

from setuptools import setup

setup(
    name='snppipeline',
    version='0.1.1',
    description='Script and functions for SNP matrix construction',
    author='Yan Luo',
    author_email='yan.luo@fda.hhs.gov',
    #url='http://', #TODO
    packages=['snppipeline'],
    long_description="""
    snppipeline is a pipeline for the production of SNP matrices from 
    sequence data used in the phylogenetic analysis of pathogenic
    organisms sequenced from samples of interest to food safety.
    """,
    #download_url='',#TODO    
#        exclude_package_data = {
#            #Exclude any ARCHIVE directories
#            '': ['ARCHIVE'],
#            #Exclude miscellaneous development notes
#            '': ['README_developmentNotes'],
#        }
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License (MIT)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords='bioinformatics NGS',
    license='GPL',
    install_requires=[
        'biopython',
        'PyVCF',
        'setuptools',
    ],
)
