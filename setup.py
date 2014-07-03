#!/usr/bin/env python2.7

from setuptools import setup
import sys

# Do not pollute innocent user's site-package with our unit tests.
# Detect the setup mode to control whether the unit test package is installed.
# Is there a better way to do this?
install_unit_tests = False
for a in sys.argv:
    if a == "test" or a == "develop":
        install_unit_tests = True

if install_unit_tests:
    packages_to_install = ['snppipeline','test']
else:
    packages_to_install = ['snppipeline']

setup(
    name='snp-pipeline',
    version='0.1.0',
    description='Script and functions for SNP matrix construction',
    author='Hugh A. Rand',
    author_email='hugh.rand@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/snp-pipeline',
    packages=packages_to_install,
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
        'scripts/prepReference.sh',
        'scripts/prepSamples.sh',
        'scripts/create_snp_matrix.py',
        'scripts/copy_snppipeline_data.py'
    ],

    # Include the test data files listed below in the distribution.
    # The package_data parameter only works for binary distributions.
    # The same list of files is in MANIFEST.in for sdist distributions.
    package_data={
        'snppipeline' : ['data/lambdaVirusInputs/reference/*.fasta',
                         'data/lambdaVirusInputs/samples/sample*/*.fastq',
                         'data/lambdaVirusExpectedResults/*.fasta',
                         'data/lambdaVirusExpectedResults/*.txt',
                         'data/lambdaVirusExpectedResults/samples/sample*/*.pileup',]
    },

    keywords=['bioinformatics', 'NGS', 'SNP'],
    license='BSD',
    install_requires=[
        'PyVCF',
        'setuptools',
    ],

    # package (aka directory) containing unit test modules
    test_suite='test',
)
