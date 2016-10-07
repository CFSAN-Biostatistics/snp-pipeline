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
will create a file called ``snppipeline.conf``::

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

SnpPipeline_StopOnSampleError
-----------------------------
Controls whether the pipeline exits upon detecting errors affecting only a single
sample.  The pipeline will always stop upon detecting global errors affecting all
samples.

**Default**: 

    When this parameter is not set to a value, the pipeline will stop upon detecting 
    single sample errors.  If you want the pipeline to continue, you must explicitly set
    this parameter false.

**Example**::

    SnpPipeline_StopOnSampleError=false


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


MaxConcurrentCallConsensus
--------------------------

Controls the number of call_consensus.py processes running concurrently 
on a workstation.  This parameter is ignored when running the pipeline on an HPC job queue.
This parameter is used by run_snp_pipeline.sh only.

**Default**: 

    When this parameter is not set to a value, the pipeline will launch multiple concurrent 
    processes using all available CPU cores on a workstation.

**Example**::

    MaxConcurrentCallConsensus=4


MaxConcurrentCollectSampleMetrics
----------------------------------

Controls the number of collectSampleMetrics.sh processes running concurrently 
on a workstation.  This parameter is ignored when running the pipeline on an HPC job queue.
This parameter is used by run_snp_pipeline.sh only.

**Default**: 

    When this parameter is not set to a value, the pipeline will launch multiple concurrent 
    processes using all available CPU cores on a workstation.

**Example**::

    MaxConcurrentCollectSampleMetrics=4


SnpPipeline_MaxSnps
-------------------
Controls the maximum number of snps allowed for each sample. Any sample with excessive snps exceeding
this limit will be excluded from the snp list, snp matrix, and snpma.vcf file. When set to -1, this 
parameter is disabled.

**Default**: 

    Do not leave this parameter unset.  To disable the excessive snp filtering and include all samples
    regardless of the number of snps, set the parameter to -1

**Example**::

    SnpPipeline_MaxSnps=1000



SnpPipeline_Aligner
-------------------
Controls which reference-based aligner is used to map reads to the reference genome.
The choices are ``bowtie2`` or ``smalt``.

**Default**: 
    
    When this parameter is not set to a value, the pipeline will use the bowtie2 aligner.

**Example**::
    
    SnpPipeline_Aligner="smalt"


Bowtie2Build_ExtraParams
------------------------

Specifies options passed to the bowtie2 indexer.  Any of the bowtie2-build options
can be specified.

**Default**: none

**Example**::

    Bowtie2Build_ExtraParams="--offrate 3"


SmaltIndex_ExtraParams
------------------------

Specifies options passed to the smalt indexer.  Any of the smalt index options
can be specified.

**Default**: none

**Example**::

    SmaltIndex_ExtraParams="-k 20 -s 1"


SamtoolsFaidx_ExtraParams
-------------------------

Specifies options passed to the SAMtools faidx indexer.  Any of the SAMtools faidx options
can be specified.

**Default**: none

**Example**::

    SamtoolsFaidx_ExtraParams=""


Bowtie2Align_ExtraParams
------------------------

Specifies options passed to the bowtie2 aligner.  Any of the bowtie2 aligner options
can be specified.

**Default**: 

|   If you do not specify the ``-p`` option, it defaults to 8 threads on an HPC or all cpu cores otherwise.
|      There is no way to completely suppress the -p option.
|   If Bowtie2Align_ExtraParams is not set to any value, the ``--reorder`` option is enabled by default.
|      Any value, even a single space, will suppress this default option.
|

**Parameter Notes**:

| ``-p``        : bowtie2 uses the specified number of parallel search threads
| ``--reorder`` : generate output records in the same order as the reads in the input file
| ``-X``        : maximum inter-mate fragment length for valid concordant paired-end alignments
|

**Example**::

    Bowtie2Align_ExtraParams="--reorder -p 16 -X 1000"


SmaltAlign_ExtraParams
----------------------

Specifies options passed to the smalt mapper.  Any of the smalt map options
can be specified.

**Default**: 

|   If you do not specify the ``-n`` option, it defaults to 8 threads on an HPC or all cpu cores otherwise.
|      There is no way to completely suppress the -n option.
|   If SmaltAlign_ExtraParams is not set to any value, the ``-O`` option is enabled by default.
|      Any value, even a single space, will suppress this default option.
|

**Parameter Notes**:

| ``-n`` : number of parallel alignment threads
| ``-O`` : generate output records in the same order as the reads in the input file
| ``-i`` : maximum insert size for paired-end reads
| ``-r`` : random number seed, if seed < 0 reads with multiple best mappings are reported as 'not mapped'
| ``-y`` : filters output alignments by a threshold in the number of exactly matching nucleotides
|

**Example**::

    SmaltAlign_ExtraParams="-O -i 1000 -r 1"


SamtoolsSamFilter_ExtraParams
-----------------------------
Specifies options passed to the SAMtools view tool when filtering the SAM file.  
Any of the SAMtools view options can be specified.

**Default**: 

| If SamtoolsSamFilter_ExtraParams is not set, the "-F 4" option is enabled by default.  
|    Any value, even a single space, will suppress the -F option.
|

**Parameter Notes**:

| ``-F 4``      : discard unmapped reads
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

| ``-q``    : minimum mapping quality for an alignment to be used
| ``-Q``    : minimum base quality for a base to be considered 
|

**Example**::

    SamtoolsMpileup_ExtraParams="-q 0 -Q 13"


VarscanMpileup2snp_ExtraParams
------------------------------
Specifies options passed to the Varscan mpileup2snp tool.
Any of the Varscan mpileup2snp options can be specified.

**Default**: None

**Parameter Notes**:

| ``--min-avg-qual`` : minimum base quality at a position to count a read
| ``--min-var-freq`` : minimum variant allele frequency threshold
|

**Example**::

    VarscanMpileup2snp_ExtraParams="--min-avg-qual 15 --min-var-freq 0.90"


VarscanJvm_ExtraParams
----------------------    
Specifies options passed to the Varscan Java Virtual Machine.  
Any of the JVM options can be specified.

**Default**: None

**Parameter Notes**:

| ``-Xmx300m``  : use 300 MB memory (modify as needed)
|

**Example**::

    VarscanJvm_ExtraParams="-Xmx300m"


RemoveAbnormalSnp_ExtraParams
------------------------------
Specifies options passed to the snp_filter.py script.

**Default**: None

**Parameter Notes**:

``--edge_length``
  The length of the edge regions in a contig, in which all SNPs will be removed.
``--window_size``
  The length of the window in which the number of SNPs should be no more than max_num_snp.
``--max_snp``
  The maximum number of SNPs allowed in a window.
``--out_group``
    Relative or absolute path to the file indicating outgroup samples, one sample ID per line.

**Example**::

    RemoveAbnormalSnp_ExtraParams="--edge_length 500 --window_size 1000 --max_snp 3 --out_group /path/to/outgroupSamples.txt"


CreateSnpList_ExtraParams
-------------------------
Specifies options passed to create_snp_list.py.

**Default**: None

**Example**::

    CreateSnpList_ExtraParams="--verbose 1"


CallConsensus_ExtraParams
-------------------------
Specifies options passed to call_consensus.py.

**Default**: None

**Parameter Notes**:

``--minBaseQual``
    Mimimum base quality score to count a read. All other snp filters take effect after the low-quality reads 
    are discarded.
``--minConsFreq``
    Consensus frequency. Mimimum fraction of high-quality reads supporting the consensus to make a call.
``--minConsStrdDpth``
    Consensus strand depth. Minimum number of high-quality reads supporting the consensus which must be present 
    on both the forward and reverse strands to make a call
``--minConsStrdBias``
    Strand bias. Minimum fraction of the high-quality consensus-supporting reads which must be present on both 
    the forward and reverse strands to make a call. The numerator of this fraction is the number of high-quality 
    consensus-supporting reads on one strand. The denominator is the total number of high-quality 
    consensus-supporting reads on both strands combined.
``--vcfFileName``
    VCF Output file name. If specified, a VCF file with this file name will be created in the same directory
    as the consensus fasta file for this sample.
``--vcfAllPos``
    Flag to cause VCF file generation at all positions, not just the snp positions. This has no effect on the
    consensus fasta file, it only affects the VCF file. This capability is intended primarily as a diagnostic
    tool and enabling this flag will greatly increase execution time.
``--vcfPreserveRefCase``
    Flag to cause the VCF file generator to emit each reference base in uppercase/lowercase as it appears in the 
    reference sequence file.  If not specified, the reference bases are emitted in uppercase.

**Example**::

    CallConsensus_ExtraParams="--verbose 1 --minBaseQual 15 --vcfFileName consensus.vcf"


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


MergeVcf_ExtraParams
--------------------
Specifies options passed to mergeVcf.sh

**Default**: none

**Example**::

    MergeVcf_ExtraParams="-n sample.vcf"


CollectSampleMetrics_ExtraParams
--------------------------------
Specifies options passed to collectSampleMetrics.sh

**Default**: none

**Example**::

    CollectSampleMetrics_ExtraParams="-v consensus.vcf"


CombineSampleMetrics_ExtraParams
--------------------------------
Specifies options passed to combineSampleMetrics.sh

**Default**: none

**Parameter Notes**:

| ``-s``  : Emit column headings with spaces instead of underscores
|

**Example**::

    CombineSampleMetrics_ExtraParams="-s"


Torque_StripJobArraySuffix
--------------------------
Controls stripping the suffix from the job id when specifying Torque job array dependencies.
It may be necessary to change this parameter if run_snp_pipeline.sh fails with an illegal qsub 
dependency error.

**Example**::

    Torque_StripJobArraySuffix=false


GridEngine_StripJobArraySuffix
------------------------------
Controls stripping the suffix from the job id when specifying Grid Engine job array dependencies.
It may be necessary to change this parameter if run_snp_pipeline.sh fails with an illegal qsub 
dependency error.

**Example**::

    GridEngine_StripJobArraySuffix=true


GridEngine_PEname
-----------------
Specifies the name of the Grid Engine parallel environment.  This is only needed when running
the SNP Pipeline on a High Performance Computing cluster with the Grid Engine job manager.  
Contact your HPC system administrator to determine the name of your parallel environment. 
Note: the name of this parameter was PEname in releases prior to 0.4.0.

**Example**::

    GridEngine_PEname="mpi"


GridEngine_QsubExtraParams
--------------------------
Specifies extra options passed to qsub when running the SNP Pipeline on the Grid Engine job scheduler.

**Default**: None

**Example**::

    GridEngine_QsubExtraParams="-q bigmem.q"


Torque_QsubExtraParams
--------------------------
Specifies extra options passed to qsub when running the SNP Pipeline on the Torque job scheduler.

**Default**: None

**Example**::

    Torque_QsubExtraParams="-l pmem=16gb"
