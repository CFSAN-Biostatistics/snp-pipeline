===============================
CFSAN SNP Pipeline
===============================

.. Image showing the PyPI version badge - links to PyPI
.. image:: https://badge.fury.io/py/snp-pipeline.png
    :target: http://badge.fury.io/py/snp-pipeline
    
.. Image showing the Travis Continuous Integration test status
.. .. image:: https://travis-ci.org/CFSAN-Biostatistics/snp-pipeline.png?branch=master
..        :target: https://travis-ci.org/CFSAN-Biostatistics/snp-pipeline

.. Image showing the PyPi download count
.. image:: https://pypip.in/d/snp-pipeline/badge.png
        :target: https://pypi.python.org/pypi/snp-pipeline


The CFSAN SNP Pipeline is a Python-based system for the production of SNP 
matrices from sequence data used in the phylogenetic analysis of pathogenic 
organisms sequenced from samples of interest to food safety.

The SNP Pipeline was developed by the United States Food 
and Drug Administration, Center for Food Safety and Applied Nutrition.

* Free software: See license below. 
* Documentation: http://snp-pipeline.rtfd.org
* Source Code: https://github.com/CFSAN-Biostatistics/snp-pipeline
* PyPI Distribution: https://pypi.python.org/pypi/snp-pipeline

Features
--------

The CFSAN SNP Pipeline uses reference-based alignments to create a matrix of
SNPs for a given set of samples. The process generally starts off by finding
a reference that is appropriate for the samples of interest, and collecting
the sample sequence data into an appropriate directory structure. The SNP
pipeline can then be used to perform the alignment of the samples to the
reference. Once the sample sequences are aligned, a list of SNP positions is
generated. The list of SNP positions is then used in combination with
alignments of the samples to the reference sequence to call SNPs. The SNP
calls are organized into a matrix containing (only) the SNP calls for all
of the sequences.

This software was developed with the objective of creating high quality
SNP matrices for sequences from closely-related pathogens, e.g., different
samples of Salmonella enteriditis from an outbreak investigation. The
focus on closely related sequences means that this code is not suited for 
the analysis of relatively distantly related organisms, where there is not
a single reference sequence appropriate for all the organisms for which an
analysis is desired.

The CFSAN SNP Pipeline is written in a combination of bash and python. The
code (including the bash scripts) is designed to be straighforward to
install. A script is provided that allows the running of the Python code
from the command line. All of the parameters used are available as command
line options. In situations where additional customization is desired, the
code is not highly complex and should be easy to modify as necessary.

Two examples of using the code are provided. These examples serve as both
unit tests, and as examples that can be modified to work on other data
sets of interest.


Citing SNP Pipeline
-------------------

To cite SNP Pipeline, please reference the SNP Pipeline GitHub repository:

    https://github.com/CFSAN-Biostatistics/snp-pipeline

and cite the associated paper:

    Pettengill, JB, Luo Y, Rand H, Chen Y, Gonzalez-Escalona N, Ottesen A,
    Allard MW, Strain E. *In review*. An evaluation of alternative methods for
    constructing phylogenies from whole genome sequence data: A case
    study with *Salmonella*. PeerJ.


License
-------

See the LICENSE.txt file included in the SNP Pipeline distribution.

