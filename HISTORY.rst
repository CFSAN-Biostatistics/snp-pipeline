.. :changelog:

History
-------

0.1.1 (2014-07-??)
++++++++++++++++++

* Moved the bowtie alignment to a new script, alignSampleToReference.sh, for 
better utilization of CPU cores when running in HPC environment.
* Changed the prepSamples.sh calling convention to take the sample directory,
not the sample files.
* prepSamples.sh uses the CLASSPATH environment variable to locate VarScan.jar.

0.1.0 (2014-07-03)
++++++++++++++++++

* Basic functionality implemented.
* Lambda virus tests created and pass.
* S. Agona tests created -- UNDER DEVELOPMENT
* Installs properly from PyPI.
* Documentation available at ReadTheDocs.
