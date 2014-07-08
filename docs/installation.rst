.. _installation-label:

============
Installation
============

.. highlight:: bash

The SNP Pipeline software package consists of python scripts and shell scripts
with dependencies on executable programs launched by the scripts.

Operating System Requirements
-----------------------------

The SNP Pipeline runs in a Linux environment. It has been tested 
on the following platforms:

    * CentOS
    * Ubuntu
    * Mac OS X

Executable Software Dependencies
--------------------------------

You should have the following software installed before using the SNP Pipeline.

    * Bowtie2_, a tool for aligning reads to long reference sequences.
    * SAMtools_, utilities for manipulating alignments in the SAM format.
    * VarScan_, a tool to detect variants in NGS data.
    * fastq-dump_, an SRA Toolkit utility for fetching samples from NCBI SRA.

Environment Variables
---------------------

Define the CLASSPATH environment variable to specify the location of the VarScan jar file.  Add 
the following (or something similiar) to your .bashrc file::

    export CLASSPATH=~/software/varscan.v2.3.6/VarScan.v2.3.6.jar:$CLASSPATH



Installing Python Package Dependencies
--------------------------------------

For the most part, the installer automatically installs the necessary python packages used by snp-pipeline.  However, not all python packages can be reliably installed automatically.  The packages listed below will need to be manually installed if not already provided by your python distribution.

    * Biopython_, a set of tools for biological computation written in Python.

Installing the SNP Pipeline Python Package
------------------------------------------
There is more than one way to install the SNP Pipeline depending on whether you intend to work with the source code or just run it.

Installation Method 1 for Software Developers
`````````````````````````````````````````````

If you intend to work with the source code in the role of a software developer, you should clone the GitHub repository as described in the :ref:`contributing-label` section of this documentation.

Installation Method 2 for Everyone Else
```````````````````````````````````````

If you merely want to run the software without viewing or changing the source code, follow the instructions below.

At the command line::

    $ pip install snp-pipeline

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv snp-pipeline
    $ pip install snp-pipeline


Uninstalling SNP Pipeline 
-------------------------

If you installed with pip, you can uninstall from the command line::

    $ pip uninstall snp-pipeline

.. _Bowtie2: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
.. _SAMtools: http://sourceforge.net/projects/samtools/files/
.. _VarScan: http://sourceforge.net/projects/varscan/files/
.. _fastq-dump: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
.. _Biopython: http://biopython.org/wiki/Download
