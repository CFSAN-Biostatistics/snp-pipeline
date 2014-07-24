.. :changelog:

History
-------

0.1.1 (2014-07-??)
~~~~~~~~~~~~~~~~~~

Bug fixes:

* The snp list, snp matrix, and referenceSNP files were incorrectly sorted by 
  position alphabetically, not numerically.
* The SNP Pipeline produced slightly different pileups each time we ran the pipeline.  
  Often we noticed two adjacent read-bases swapped in the pileup files.  This was 
  caused by utilizing multiple CPU cores during the bowtie alignment.  The output 
  records in the SAM file were written in non-deterministic order when bowtie ran 
  with multiple concurrent threads.  Fixed by adding the ``--reorder`` option to the 
  bowtie alignment command line.

Other Changes: 

*Note the loss of backward compatibilty for existing workflows using prepSamples.sh*

* Moved the bowtie alignment to a new script, alignSampleToReference.sh, for 
  better control of CPU core utilization when running in HPC environment.
* Changed the prepSamples.sh calling convention to take the sample directory, 
  not the sample files. 
* prepSamples.sh uses the CLASSPATH environment variable to locate VarScan.jar.

0.1.0 (2014-07-03)
~~~~~~~~~~~~~~~~~~

* Basic functionality implemented.
* Lambda virus tests created and pass.
* S. Agona tests created -- UNDER DEVELOPMENT
* Installs properly from PyPI.
* Documentation available at ReadTheDocs.
