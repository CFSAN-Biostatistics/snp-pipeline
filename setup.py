#!/usr/bin/env python2.7

from setuptools import setup

setup(
    name='snp-pipeline',
    version='0.1.1',
    description='Script and functions for SNP matrix construction',
    author='Yan Luo',
    author_email='yan.luo@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/snp-pipeline',
    packages=['snppipeline','test'],
    long_description="""
    snp-pipeline is a pipeline for the production of SNP matrices from
    sequence data used in the phylogenetic analysis of pathogenic
    organisms sequenced from samples of interest to food safety.
    """,

    # download_url should be used if the distribution is not hosted on PyPI.
    # GitHub can host the downloadable distribution if the repo is tagged
    #download_url='',#TODO    #TODO set this up properly

    #TODO figure out how to next bit up
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
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    scripts=[
        'scripts/prepSamples.sh',
        'scripts/prepReference.sh',
        'scripts/prepSequenceData.sh',
        'scripts/runsnppipeline.py'
    ],

    # Include the test data files listed below in the distribution.
    # The package_data parameter only works for binary distributions.
    # The same list of files is in MANIFEST.in for sdist distributions.
    package_data={
        'test' : ['testLambdaVirus/sample*/*.fq',
                  'codeComparisonFiles/testLambdaVirus/*.fasta',
                  'codeComparisonFiles/testLambdaVirus/sample*/*.pileup']
    },

    keywords=['bioinformatics', 'NGS', 'SNP'],
    license='BSD',
    install_requires=[
        'biopython',
        'PyVCF',
        'setuptools',
    ],

    # package (aka directory) containing unit test modules
    test_suite='test',
)
