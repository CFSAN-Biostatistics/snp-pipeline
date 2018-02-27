.. _installation-label:

============
Installation
============

.. highlight:: bash

The SNP Pipeline software package consists of python scripts and shell scripts
with dependencies on executable programs launched by the scripts.

Step 1 - Operating System Requirements
--------------------------------------
The SNP Pipeline runs in a Linux environment. It has been tested
on the following platforms:

    * Red Hat
    * CentOS
    * Ubuntu

Step 2 - Executable Software Dependencies
-----------------------------------------
You should have the following software installed before using the SNP Pipeline.

    * Bowtie2_, a tool for aligning reads to long reference sequences
    * SMALT_, a tool for aligning reads to long reference sequences
    * SAMtools_, utilities for manipulating alignments in the SAM format
    * Picard_, a set of tools for manipulating sequencing data in SAM/BAM/CRAM/VCF format
    * GATK_, variant discovery and genotyping tools
    * VarScan_, a tool to detect variants in NGS data
    * tabix_, a generic indexer for tab-delimited genome position files
    * bgzip, part of the tabix package, bgzip is a block compression utility
    * BcfTools_, utilities for variant calling and manipulating VCFs and BCFs
    * fastq-dump_, an SRA Toolkit utility for fetching samples from NCBI SRA

Note: you will need either Bowtie2 or SMALT.  You do not have to install both.
However, the included result files were generated with Bowtie2.  Your results may differ
when using SMALT.

Step 3 - Environment Variables
------------------------------
Define the CLASSPATH environment variable to specify the location of the Picard and VarScan jar files.  Add
the following lines (or something similiar) to your .bashrc file::

    export CLASSPATH=~/software/varscan.v2.3.9/VarScan.jar:$CLASSPATH
    export CLASSPATH=~/software/picard/picard.jar:$CLASSPATH
    export CLASSPATH=~/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar:$CLASSPATH

The order is important.  The GATK path should be before Picard in the CLASSPATH environment variable.

Step 4 - Python
---------------
The SNP pipeline is compatible with python version 2.6, 2.7, 3.3, 3.4, and 3.5.  The pipeline has not been tested on other python versions.
If you do not already have python installed, you should install version 2.7.  You can either build from source
or install a precompiled version with your Linux package manager.


Step 5 - Pip
------------
This can be a troublesome installation step -- proceed with caution.  The pip tool is used to install python packages
including the snp-pipeline and other packages used by the snp-pipeline.  Some newer versions of Python include pip.
Check to see if pip is already installed::

    $ pip -V

If pip is not already installed, proceed as follows::

    Download get-pip.py from https://pip.pypa.io/en/latest/installing.html#install-pip
    $ python get-pip.py --user

Note: avoid using sudo when installing pip.  Some users have experienced problems installing and loading packages when pip is installed using sudo.


Step 6 - Python Package Dependencies
------------------------------------

For the most part, the installer automatically installs the necessary python packages used by snp-pipeline.  However,
not all python packages can be reliably installed automatically.  The packages listed below may need to be manually
installed if automatic installation fails.  You can either install these packages
now, or hope for the best and manually install later if the automatic installation fails.

    * Biopython_, a set of tools for biological computation written in Python.

Step 7 - Install the SNP Pipeline Python Package
------------------------------------------------
There is more than one way to install the SNP Pipeline depending on whether you intend to work with the source code or just run it.

Installation Method 1 for Most Users
````````````````````````````````````

This is the recommended installation method for new users.

If you want to run the software without viewing or changing the source code, follow the instructions below.

At the command line::

    $ pip install --user snp-pipeline

Update your .bashrc file with the path to user-installed python packages::

    export PATH=~/.local/bin:$PATH

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv snp-pipeline
    $ pip install snp-pipeline



Installation Method 2 for Software Developers
`````````````````````````````````````````````

If you intend to work with the source code in the role of a software developer, you should clone the GitHub repository as described in the :ref:`contributing-label` section of this documentation.


Upgrading SNP Pipeline
----------------------
If you previously installed with pip, you can upgrade to the newest version from the command line::

    $ pip install --user --upgrade snp-pipeline


Uninstalling SNP Pipeline
-------------------------

If you installed with pip, you can uninstall from the command line::

    $ pip uninstall snp-pipeline

Tips
----

There is a dependency on the python psutil package.  Pip will attempt to
install the psutil package automatically when installing snp-pipeline.
If it fails with an error message about missing Python.h, you will need to
manually install the python-dev package.
In Ubuntu, use this command::

    $ sudo apt-get install python-dev

You may need to upgrade your Java Runtime Environment (JRE) to run Picard.


.. _Bowtie2: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
.. _SAMtools: http://sourceforge.net/projects/samtools/files/
.. _Picard: https://broadinstitute.github.io/picard/command-line-overview.html
.. _GATK: https://software.broadinstitute.org/gatk/download/archive
.. _VarScan: http://sourceforge.net/projects/varscan/files/
.. _tabix: http://www.htslib.org/doc/tabix.html
.. _BcfTools: http://sourceforge.net/projects/samtools/files/samtools/1.1/
.. _fastq-dump: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
.. _Biopython: http://biopython.org/wiki/Download
.. _SMALT: http://sourceforge.net/projects/smalt/files
