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

TODO: Fix this

You should have the following software installed before using the SNP Pipeline.

    * Bowtie2_, a tool for aligning reads to long reference sequences.
    * SAMtools_, utilities for manipulating alignments in the SAM format.
    * BCFtools_, utilities for variant calling and manipulating VCF/BCF files.
    * fastq-dump_, an SRA Toolkit utility for fetching samples from NCBI SRA


Python Script
-------------

TODO: Fix this

At the command line::

    $ easy_install snp-pipeline

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv snp-pipeline
    $ pip install snp-pipeline


.. _Bowtie2: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
.. _SAMtools: http://sourceforge.net/projects/samtools/files/
.. _BCFtools: http://samtools.github.io/bcftools/
.. _SRA Toolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software