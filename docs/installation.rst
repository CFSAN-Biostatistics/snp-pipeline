============
Installation
============

The SNP Pipeline software package consists of python scripts and shell scripts
with dependencies on executable programs launched by the scripts.

Operating System Requirements
----------------------------

The SNP Pipeline runs in a Linux environment. It has been tested 
on the following platforms:
    * CentOS
    * Ubuntu
    * Mac OS X

Executable Software Dependencies
--------------------------------

TODO: what other software is needed?

You should have the following software installed before using the SNP Pipeline.

    * Bowtie2_, a tool for aligning reads to long reference sequences.
    * SAMtools_, utilities for manipulating alignments in the SAM format.
    * BCFtools_, utilities for variant calling and manipulating VCF/BCF files.
    * fastq-dump_, an SRA Toolkit utility for fetching samples from NCBI SRA


Installing Python Package Dependencies
--------------------------------------

The installer automatically installs the necessary python packages used by
snp-pipeline.  However, not all python packages can be reliably installed
automatically.  The packages listed below may need manual installation if
not already provided by your python distribution.

    * Biopython

Installing the SNP Pipeline Python Package
------------------------------------------
There is more than one way to install the SNP Pipeline.  If you intend to 
work with the source code in the role of a software developer, you should
clone the GitHub repository as described here:

.. include:: ../CONTRIBUTING.rst

Non-software developers should follow the instructions below.

TODO: Fix this

At the command line::

    $ easy_install snp-pipeline

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv snp-pipeline
    $ pip install snp-pipeline


.. _Bowtie2: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
.. _SAMtools: http://sourceforge.net/projects/samtools/files/
.. _BCFtools: http://samtools.github.io/bcftools/
.. _fastq-dump: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software