.. _configuration-label:

=============
Configuration
=============

You can customize the behavior of the SNP Pipeline by configuring parameters.  
Each step in the pipeline has a corresponding parameter allowing you to set one 
or more options for each tool in the pipeline. 

Parameters can be configured either in a configuration file if you are using the
``run_snp_pipeline.sh`` script, or in environment variables.

The pipeline comes with a default configuration file.  When you run the run_snp_pipeline.sh 
script without specifying a configuration file, it automatically uses the 
default supplied configuration file.  

To get a copy of the default configuration file, run the following command.  This 
will create a file called snppipeline.conf::

    copy_snppipeline_data.py configurationFile

To customize the pipeline behavior, edit the configuration file and pass the file to
the run_snp_pipeline.sh script::

    run_snp_pipeline.sh -c snppipeline.conf ...

When the run_snp_pipeline.sh script runs, it copies the configuration file to the
log directory for the run, capturing the configuration used for each run.

If you decide not to use the run_snp_pipeline.sh script, you can still customize the
behavior of the pipeline tools, but you will need to set (and export) environment 
variables with the same names as the parameters in the configuration file.

The available configuration parameters are described below.


MaxConcurrentPrepSamples
------------------------

Controls the number of prepSamples.sh (SAMtools and Varscan) processes running concurrently 
on a workstation.  This parameter is ignored when running the pipeline on an HPC job queue.
This parameter is used by run_snp_pipeline.sh only.

**Default**: 

    When this parameter is not set to a value, the pipeline will launch multiple concurrent 
    processes using all available CPU cores on a workstation.

**Example**::

    MaxConcurrentPrepSamples=2


MaxConcurrentCreateSnpPileup
----------------------------

Controls the number of create_snp_pileup.py  processes running concurrently 
on a workstation.  This parameter is ignored when running the pipeline on an HPC job queue.
This parameter is used by run_snp_pipeline.sh only.

**Default**: 

    When this parameter is not set to a value, the pipeline will launch multiple concurrent 
    processes using all available CPU cores on a workstation.

**Example**::

    MaxConcurrentCreateSnpPileup=4


MaxConcurrentCollectSampleMetrics
----------------------------------

Controls the number of collectSampleMetrics.sh  processes running concurrently 
on a workstation.  This parameter is ignored when running the pipeline on an HPC job queue.
This parameter is used by run_snp_pipeline.sh only.

**Default**: 

    When this parameter is not set to a value, the pipeline will launch multiple concurrent 
    processes using all available CPU cores on a workstation.

**Example**::

    MaxConcurrentCollectSampleMetrics=4


Bowtie2Build_ExtraParams
------------------------

Specifies options passed to the bowtie2 indexer.  Any of the bowtie2-build options
can be specified.

**Default**: none

**Example**::

    Bowtie2Build_ExtraParams="--offrate 3"


SamtoolsFaidx_ExtraParams
-------------------------

Specifies options passed to the SAMtools faidx indexer.  Any of the SAMtools faidx options
can be specified.

**Default**: none

**Example**::

    SamtoolsFaidx_ExtraParams=""


Bowtie2Align_ExtraParams
------------------------

Specifies options passed to the bowtie2 aligner indexer.  Any of the bowtie2 aligner options
can be specified.

**Default**: 

|   If you do not specify the -p option, it defaults to 8 threads on an HPC or all cpu cores otherwise.
|      There is no way to completely suppress the -p option.
|   If Bowtie2Align_ExtraParams is not set to any value, the ``--reorder`` option is enabled by default.
|      Any value, even a single space, will suppress this default option.
|

**Parameter Notes**:

| -p        : bowtie2 uses the specified number of parallel search threads
| --reorder : generate output records in the same order as the reads in the input file
|

**Example**::

    Bowtie2Align_ExtraParams="--reorder -p 16"


SamtoolsSamFilter_ExtraParams
-----------------------------
Specifies options passed to the SAMtools view tool when filtering the SAM file.  
Any of the SAMtools view options can be specified.

**Default**: 

| If SamtoolsSamFilter_ExtraParams is not set, the "-F 4" option is enabled by default.  
|    Any value, even a single space, will suppress the -F option.
|

**Parameter Notes**:

| -F 4      : discard unmapped reads
|

**Example**::

    SamtoolsSamFilter_ExtraParams="-F 4"


SamtoolsSort_ExtraParams
------------------------
Specifies options passed to the SAMtools sort tool when sorting the BAM file.  
Any of the SAMtools sort options can be specified.

**Default**: None

**Example**::

    SamtoolsSort_ExtraParams=""


SamtoolsMpileup_ExtraParams
---------------------------
Specifies options passed to the SAMtools mpileup tool.  
Any of the SAMtools mpileup options can be specified.

**Default**: None

**Parameter Notes**:

| -q    : minimum mapping quality for an alignment to be used
| -Q    : minimum base quality for a base to be considered 
|

**Example**::

    SamtoolsMpileup_ExtraParams="-q 0 -Q 13"


VarscanMpileup2snp_ExtraParams
------------------------------
Specifies options passed to the Varscan mpileup2snp tool.
Any of the Varscan mpileup2snp options can be specified.

**Default**: None

**Parameter Notes**:

| --min-avg-qual : minimum base quality at a position to count a read
| --min-var-freq : minimum variant allele frequency threshold
|

**Example**::

    VarscanMpileup2snp_ExtraParams="--min-avg-qual 15 --min-var-freq 0.90"


VarscanJvm_ExtraParams
----------------------    
Specifies options passed to the Varscan Java Virtual Machine.  
Any of the JVM options can be specified.

**Default**: None

**Parameter Notes**:

| -Xmx300m  : use 300 MB memory (modify as needed)
|

**Example**::

    VarscanJvm_ExtraParams="-Xmx300m"


CreateSnpList_ExtraParams
-------------------------
Specifies options passed to create_snp_list.py.

**Default**: None

**Example**::

    CreateSnpList_ExtraParams="--verbose 1"


CreateSnpPileup_ExtraParams
---------------------------
Specifies options passed to create_snp_pileup.py.

**Default**: None

**Example**::

    CreateSnpPileup_ExtraParams="--verbose 1"


CreateSnpMatrix_ExtraParams
---------------------------
Specifies options passed to create_snp_matrix.py.

**Default**: None

**Example**::

    CreateSnpMatrix_ExtraParams="--verbose 1"


CreateSnpReferenceSeq_ExtraParams
---------------------------------
Specifies options passed to create_snp_reference_seq.py.

**Default**: None

**Example**::

    CreateSnpReferenceSeq_ExtraParams="--verbose 1"

PEname
------
Specifies the name of the Grid Engine parallel environment.  This is only needed when running
the SNP Pipeline on a High Performance Computing cluster with the Grid Engine job manager.  
Contact your HPC system administrator to determine the name of your parallel environment.

**Example**::

    PEname="mpi"
