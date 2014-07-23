.. :changelog:

History
-------

0.1.1 (2014-07-??)
~~~~~~~~~~~~~~~~~~

Bug fixes:

* The snp list, snp matrix, and referenceSNP files were incorrectly sorted by 
  position alphabetically, not numerically.

Other Changes:

* Moved the bowtie alignment to a new script, alignSampleToReference.sh, for 
  better control of CPU core utilization when running in HPC environment.
* Changed the prepSamples.sh calling convention to take the sample directory, 
  not the sample files. Note: this breaks backward compatibilty for any scripts
  written with v0.1.0 prepSamples.sh
* prepSamples.sh uses the CLASSPATH environment variable to locate VarScan.jar.

0.1.0 (2014-07-03)
~~~~~~~~~~~~~~~~~~

* Basic functionality implemented.
* Lambda virus tests created and pass.
* S. Agona tests created -- UNDER DEVELOPMENT
* Installs properly from PyPI.
* Documentation available at ReadTheDocs.
