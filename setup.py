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


# Control which 3rd party packages should be installed 
# depending on the python version
install_requires = [
    'PyVCF',
    'setuptools',
    'psutil',
    'Biopython',
]

# Below needed for Python 2.6
if sys.version_info < (2,7,):
    install_requires.append('argparse')
    install_requires.append('ordereddict')
    install_requires.append('counter')

test_requires = [
    'testfixtures',
]

setup(
    name='snp-pipeline',
    version='0.4.0',
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
        'scripts/run_snp_pipeline.sh',
        'scripts/prepReference.sh',
        'scripts/alignSampleToReference.sh',
        'scripts/prepSamples.sh',
        'scripts/create_snp_list.py',
        'scripts/create_snp_pileup.py',
        'scripts/call_consensus.py',
        'scripts/create_snp_matrix.py',
        'scripts/create_snp_reference_seq.py',
        'scripts/copy_snppipeline_data.py',
        'scripts/collectSampleMetrics.sh',
        'scripts/combineSampleMetrics.sh',
        'scripts/mergeVcf.sh',
    ],

    # Include the test data files listed below in the distribution.
    # The package_data parameter only works for binary distributions.
    # The same list of files is in MANIFEST.in for sdist distributions.
    package_data={
        'snppipeline' : ['data/configuration/snppipeline.conf',
                         'data/lambdaVirusInputs/reference/*',
                         'data/lambdaVirusInputs/samples/*/*',
                         'data/lambdaVirusExpectedResults/samples/*/*',
                         'data/lambdaVirusExpectedResults/*.*',
                         'data/agonaInputs/reference/*',
                         'data/agonaInputs/sha256sumCheck',
                         'data/agonaExpectedResults/*.*',
                         'data/agonaExpectedResults/samples/*/*',
                         'data/listeriaInputs/reference/*',
                         'data/listeriaInputs/sampleList',
                         'data/listeriaInputs/sha256sumCheck',
                         'data/listeriaExpectedResults/*.*',
                         'data/listeriaExpectedResults/samples/*/*',
                         ]
    },

    keywords=['bioinformatics', 'NGS', 'SNP'],
    license='BSD',
    install_requires=install_requires,

    # package (aka directory) containing unit test modules
    test_suite='test',
    tests_require=test_requires,
)
