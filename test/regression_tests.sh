#!/bin/bash
#
# Author: Steven C. Davis (scd)
# Purpose: SNP Pipline regression tests
# Input:
#       none
# Output:
#       Test results written to stdout
# Use example:
#       ./regression_tests.sh
# History:
#       20150417-scd: Started.
# Notes:
#
# Bugs:
#
# References:
#       https://code.google.com/p/shunit2/wiki/ProjectInfo
#       http://shunit2.googlecode.com/svn/trunk/source/2.1/doc/shunit2.html
#

#################################################
# Helper Functions
#################################################

# Fail the current test with an error message if a file is missing or not empty
verifyEmptyFile()
{
    if [[ ! -f "$1" ]]; then fail "The file $1 does not exist."; return $SHUNIT_FALSE; fi
    if [[ -s "$1" ]]; then fail "The file $1 is not empty.";  return $SHUNIT_FALSE; fi
}

# Fail the current test with an error message if a file is missing, empty, or not readable.
verifyNonEmptyReadableFile()
{
    if [[ ! -f "$1" ]]; then fail "The file $1 does not exist."; return $SHUNIT_FALSE; fi
    if [[ ! -s "$1" ]]; then fail "The file $1 is empty.";  return $SHUNIT_FALSE; fi
    if [[ ! -r "$1" ]]; then fail "The file $1 is not readable.";  return $SHUNIT_FALSE; fi
}

# Fail the current test if a specified file exists
verifyNonExistingFile()
{
    if [[ -f "$1" ]]; then fail "The file $1 should not exist."; return $SHUNIT_FALSE; fi
}

# Assert a string is contained in a file.
assertFileContains()
{
    file=$1
    targetString="$2"
    grepResult=$(grep --max-count=1 -E "$targetString" "$file")
    assertNotNull "$targetString is missing from $file" "$grepResult"
}

# Assert a string is not contained in a file.
assertFileNotContains()
{
    file=$1
    targetString="$2"
    grepResult=$(grep "$targetString" "$file")
    assertNull "$targetString should not be found in $file" "$grepResult"
}

# Assert two files are identical.
# Options 3 and 4 are used to exclude some lines from the comparison with --ignore-matching-lines.
assertIdenticalFiles()
{
    if [[ ! -f "$2" ]]; then fail "The file $2 does not exist."; return $SHUNIT_FALSE; fi
    diff -q $3 $4 $5 "$1" "$2" &> /dev/null
    assertTrue "$1 does not exactly match $2" $?
}

# Assert two files are different.
# Options 3 and 4 are used to exclude some lines from the comparison with --ignore-matching-lines.
assertDifferentFiles()
{
    if [[ ! -f "$2" ]]; then fail "The file $2 does not exist."; return $SHUNIT_FALSE; fi
    diff -q $3 $4 $5 "$1" "$2" &> /dev/null
    assertFalse "The file $1 should not match $2" $?
}

# Assert file $1 is newer than file $2.
assertNewerFile()
{
    [ ! "$1" -ot "$2" ]  # $1 not older than $2
    errorCode=$?
    assertTrue "The file $1 should be newer than $2." $errorCode
}


#################################################
# Tests
#################################################


# Verify the scripts were properly installed on the path.
testScriptsOnPath()
{
    assertNotNull "cfsan_snp_pipeline is not on the path"           "$(which cfsan_snp_pipeline)"
    assertNotNull "copy_snppipeline_data.py is not on the path"     "$(which copy_snppipeline_data.py)"
    assertNotNull "run_snp_pipeline.sh is not on the path"          "$(which run_snp_pipeline.sh)"
    assertNotNull "prepReference.sh is not on the path"             "$(which prepReference.sh)"
    assertNotNull "alignSampleToReference.sh is not on the path"    "$(which alignSampleToReference.sh)"
    assertNotNull "prepSamples.sh is not on the path"               "$(which prepSamples.sh)"
    assertNotNull "snp_filter.py is not on the path"                "$(which snp_filter.py)"
    assertNotNull "create_snp_list.py is not on the path"           "$(which create_snp_list.py)"
    assertNotNull "call_consensus.py is not on the path"            "$(which call_consensus.py)"
    assertNotNull "create_snp_matrix.py is not on the path"         "$(which create_snp_matrix.py)"
    assertNotNull "create_snp_reference_seq.py is not on the path"  "$(which create_snp_reference_seq.py)"
    assertNotNull "collectSampleMetrics.sh is not on the path"      "$(which collectSampleMetrics.sh)"
    assertNotNull "combineSampleMetrics.sh is not on the path"      "$(which combineSampleMetrics.sh)"
    assertNotNull "calculate_snp_distances.py is not on the path"   "$(which calculate_snp_distances.py)"
}

# Verify cfsan_snp_pipeline data emits lambda test data
testCopySnpPipelineLambdaData()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    cfsan_snp_pipeline data lambdaVirusInputs $tempDir
    verifyNonEmptyReadableFile "$tempDir/reference/lambda_virus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/sample1_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/sample1_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/sample2_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/sample2_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/sample3_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/sample3_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/sample4_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/sample4_2.fastq"

    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix.tsv"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise_preserved.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix_preserved.tsv"
}

# Verify cfsan_snp_pipeline data emits configuration file
testCopySnpPipelineConfigurationFile()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    cfsan_snp_pipeline data configurationFile $tempDir
    verifyNonEmptyReadableFile "$tempDir/snppipeline.conf"
}


# Helper function used to check if the snp pipeline complains when it should about tools not on path
verifyWhetherCommandOnPathChecked()
{
    file=$1
    targetCommand="$2"
    which "$targetCommand" > /dev/null
    if [[ $? != 0 ]]; then
        assertFileContains "$file" "$targetCommand is not on the path"
    fi
}



# Verify run_snp_pipeline checks for necessary scripts and tools
tryRunSnpPipelineDependencyRaiseFatalError()
{
    expectErrorCode=$1
    aligner=$2

    # Save path and classpath then destroy them
    savePath="$PATH"
    saveClassPath="$CLASSPATH"
    export PATH="$(dirname $(which which)):$(dirname $(which mkdir)):$(dirname $(which run_snp_pipeline.sh))"
    export CLASSPATH=""

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Restore paths
    PATH="$savePath"
    CLASSPATH="$saveClassPath"

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the dependencies were not setup properly." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "cfsan_snp_pipeline"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "$aligner"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "samtools"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "java"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "tabix"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bgzip"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bcftools"
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to VarScan"
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to Picard"
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to GATK"
    assertFileContains "$tempDir/error.log" "Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html"

    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "cfsan_snp_pipeline"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "$aligner"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "samtools"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "java"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "tabix"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "bgzip"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "bcftools"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to VarScan"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to Picard"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to GATK"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html"
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify the index_ref command detects a misconfigured environment variable
tryIndexRefEnvironmentRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately misconfigure the environment
    export SnpPipeline_Aligner=garbage

    # Run cfsan_snp_pipeline index_ref
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref failed"
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/indexRef.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEnvironmentRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryIndexRefEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEnvironmentRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryIndexRefEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEnvironmentRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryIndexRefEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
tryIndexRefEmptyFastaFileRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately create empty fasta file
    rm "$tempDir/reference/lambda_virus.fasta"
    touch "$tempDir/reference/lambda_virus.fasta"

    # Run cfsan_snp_pipeline index_ref
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the reference file was empty." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref failed"
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileContains "$logDir/indexRef.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEmptyFastaFileRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryIndexRefEmptyFastaFileRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEmptyFastaFileRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryIndexRefEmptyFastaFileRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testIndexRefEmptyFastaFileRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryIndexRefEmptyFastaFileRaiseGlobalError 100
}


# Verify the index_ref command detects bowtie error and emits the global error marker file.
tryIndexRefBowtieIndexTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"

    # Run cfsan_snp_pipeline index_ref
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref / bowtie returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/error.log" "error|failed"
    assertFileNotContains "$logDir/indexRef.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testIndexRefBowtieIndexTrapStop()
{
    export StopOnSampleError=true
    tryIndexRefBowtieIndexTrap 100
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testIndexRefBowtieIndexTrapNoStop()
{
    export StopOnSampleError=false
    tryIndexRefBowtieIndexTrap 100
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testIndexRefBowtieIndexTrapStopUnset()
{
    unset StopOnSampleError
    tryIndexRefBowtieIndexTrap 100
}


# Verify the index_ref command detects smalt error and emits the global error marker file.
tryIndexRefSmaltIndexTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"
    sed -i 's/A/>\n/g' "$tempDir/reference/lambda_virus.fasta"

    # Run cfsan_snp_pipeline index_ref
    export SnpPipeline_Aligner=smalt
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref / smalt returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/indexRef.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/error.log" "smalt index"
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testIndexRefSmaltIndexTrapStop()
{
    export StopOnSampleError=true
    tryIndexRefSmaltIndexTrap 100
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testIndexRefSmaltIndexTrapNoStop()
{
    export StopOnSampleError=false
    tryIndexRefSmaltIndexTrap 100
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testIndexRefSmaltIndexTrapStopUnset()
{
    unset StopOnSampleError
    tryIndexRefSmaltIndexTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
tryIndexRefSamtoolsFaidxTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file, changing line lengths
    sed -i 's/A/AA/g' "$tempDir/reference/lambda_virus.fasta"

    # Run cfsan_snp_pipeline index_ref
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    # Verify error handling behavior
    #
    # The results are different depending on which version of SAMtools is installed because SAMtools 1.3
    # does not fail with an error code that can be trapped.  The verification tests below check for the
    # common results in v0.1.19 and v1.3.

    # SAMtools 0.1.19 error.log
    # -------------------------
    # Error detected while running cfsan_snp_pipeline index_ref.
    #
    # The command line was:
    #     cfsan_snp_pipeline index_ref /tmp/shunit.dlWnpz/tmp/tmp.KW30guVt3J/reference/lambda_virus.fasta
    #
    # The command at line 173 returned error code 139:
    #    samtools faidx $SamtoolsFaidx_ExtraParams "$referenceFilePath"

    # SAMtools 1.3 error.log
    # ----------------------
    # cfsan_snp_pipeline index_ref failed.
    # Error: /tmp/shunit.VB9zhq/tmp/tmp.qQsT4ORbSy/reference/lambda_virus.fasta.fai does not exist after running samtools faidx.

    # SAMtools 0.1.19 indexRef.log
    # ---------------------------------
    # [fai_build_core] different line length in sequence 'gi|9626243|ref|NC_001416.1|'.
    # /home/steven.davis/.virtualenvs/snp-pipeline-2.7/bin/prepReference.sh: line 176: 30719 Segmentation fault      (core dumped) samtools faidx $SamtoolsFaidx_ExtraParams "$referenceFilePath"

    # SAMtools 1.3 indexRef.log
    # --------------------------
    # [fai_build_core] different line length in sequence 'gi|9626243|ref|NC_001416.1|'.
    # Error: /tmp/shunit.VB9zhq/tmp/tmp.qQsT4ORbSy/reference/lambda_virus.fasta.fai does not exist after running samtools faidx.

    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/error.log" "samtools faidx"
    assertFileNotContains "$logDir/indexRef.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testIndexRefSamtoolsFaidxTrapStop()
{
    export StopOnSampleError=true
    tryIndexRefSamtoolsFaidxTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testIndexRefSamtoolsFaidxTrapNoStop()
{
    export StopOnSampleError=false
    tryIndexRefSamtoolsFaidxTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testIndexRefSamtoolsFaidxTrapUnset()
{
    unset StopOnSampleError
    tryIndexRefSamtoolsFaidxTrap 100
}


# Verify the index_ref command detects picard CreateSequenceDictionary error and emits the global error marker file.
tryIndexRefCreateSequenceDictionaryTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables cfsan_snp_pipeline index_ref will use.
    # This simulates what run_snp_pipeline does before running cfsan_snp_pipeline index_ref.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline index_ref, so all prior steps will be successful
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"

    # Remove the dict file, so it will be created
    rm "$tempDir/reference/lambda_virus.dict"

    # Deliberately corrupt the fasta and fai file -- CreateSequenceDictionary will use fai file if it exists
    head -c 10000 /dev/urandom > "$tempDir/reference/lambda_virus.fasta"
    head -c 10000 /dev/urandom > "$tempDir/reference/lambda_virus.fasta.fai"

    # Adjust file timestamps to skip bowtie and samtools faidx
    touch -d  '-3 day' "$tempDir/reference/lambda_virus.fasta"
    touch -d  '-2 day' "$tempDir/reference/lambda_virus*bt2"
    touch -d  '-1 day' "$tempDir/reference/lambda_virus*fai"

    # Run cfsan_snp_pipeline index_ref again
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    errorCode=$?

    assertEquals "cfsan_snp_pipeline index_ref CreateSequenceDictionary returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/error.log" "CreateSequenceDictionary"
    assertFileNotContains "$logDir/indexRef.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/indexRef.log" "lambda_virus.dict is already freshly built"
}

# Verify the index_ref command detects picard CreateSequenceDictionary error and emits the global error marker file.
testIndexRefCreateSequenceDictionaryTrapStop()
{
    export StopOnSampleError=true
    tryIndexRefCreateSequenceDictionaryTrap 100
}

# Verify the index_ref command detects picard CreateSequenceDictionary error and emits the global error marker file.
testIndexRefCreateSequenceDictionaryTrapNoStop()
{
    export StopOnSampleError=false
    tryIndexRefCreateSequenceDictionaryTrap 100
}

# Verify the index_ref command detects picard CreateSequenceDictionary error and emits the global error marker file.
testIndexRefCreateSequenceDictionaryTrapUnset()
{
    unset StopOnSampleError
    tryIndexRefCreateSequenceDictionaryTrap 100
}


# Verify the map_reads command detects a misconfigured environment variable.
tryMapReadsEnvironmentRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"

    # Deliberately misconfigure the environment
    export SnpPipeline_Aligner=garbage

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/mapReads.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the map_reads command detects a misconfigured environment variable.
testMapReadsEnvironmentRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsEnvironmentRaiseGlobalError 100
}

# Verify the map_reads command detects a misconfigured environment variable.
testMapReadsEnvironmentRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsEnvironmentRaiseGlobalError 100
}

# Verify the map_reads command detects a misconfigured environment variable.
testMapReadsEnvironmentRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsEnvironmentRaiseGlobalError 100
}


# Verify the map_reads command detects an Missing reference file
tryMapReadsMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/mapReads.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Use the -f option to force a rebuild"
}

# Verify the map_reads command detects a missing reference file
testMapReadsMissingReferenceRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsMissingReferenceRaiseGlobalError 100
}

# Verify the map_reads command detects a missing reference file
testMapReadsMissingReferenceRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsMissingReferenceRaiseGlobalError 100
}

# Verify the map_reads command detects a missing reference file
testMapReadsMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsMissingReferenceRaiseGlobalError 100
}


# Verify the map_reads command detects a missing sample file
tryMapReadsMissingSample1RaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileContains "$logDir/mapReads.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Use the -f option to force a rebuild"
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testMapReadsMissingSample1RaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsMissingSample1RaiseSampleError 100
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testMapReadsMissingSample1RaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsMissingSample1RaiseSampleError 98
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testMapReadsMissingSample1RaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsMissingSample1RaiseSampleError 100
}


# Verify the map_reads command detects a missing sample file
tryMapReadsMissingSample2RaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/samples/sample1/sample1_2.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileContains "$logDir/mapReads.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testMapReadsMissingSample2RaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsMissingSample2RaiseSampleError 100
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testMapReadsMissingSample2RaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsMissingSample2RaiseSampleError 98
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testMapReadsMissingSample2RaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsMissingSample2RaiseSampleError 100
}


# Verify the map_reads command detects bowtie alignment error.
tryMapReadsBowtieAlignTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"
    export SnpPipeline_Aligner=bowtie2

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log-1"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log-1" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/mapReads.log-1" "bowtie2"
    assertFileNotContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/mapReads.log-3"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log-3" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/mapReads.log-3" "bowtie2"
    assertFileNotContains "$logDir/mapReads.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log-3" "Use the -f option to force a rebuild"
}


# Verify the map_reads command detects bowtie alignment error.
testMapReadsBowtieAlignTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsBowtieAlignTrap 100
}

# Verify the map_reads command detects bowtie alignment error.
testMapReadsBowtieAlignTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsBowtieAlignTrap 98
}

# Verify the map_reads command detects bowtie alignment error.
testMapReadsBowtieAlignTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsBowtieAlignTrap 100
}


# Verify the map_reads command detects smalt alignment error.
tryMapReadsSmaltAlignTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"
    export SnpPipeline_Aligner=smalt

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log-1"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log-1" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/mapReads.log-1" "smalt map"
    assertFileNotContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/mapReads.log-3"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log-3" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/mapReads.log-3" "smalt map"
    assertFileNotContains "$logDir/mapReads.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log-3" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}


# Verify the map_reads command detects smalt alignment error.
testMapReadsSmaltAlignTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsSmaltAlignTrap 100
}

# Verify the map_reads command detects smalt alignment error.
testMapReadsSmaltAlignTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsSmaltAlignTrap 98
}

# Verify the map_reads command detects smalt alignment error.
testMapReadsSmaltAlignTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsSmaltAlignTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects a Missing reference file
tryCallSitesMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try cfsan_snp_pipeline call_sites with missing file
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta" "xxxx" &> "$logDir/callSites.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/callSites.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/callSites.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing Reference file
testCallSitesMissingReferenceRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryCallSitesMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing reference file
testCallSitesMissingReferenceRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryCallSitesMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing reference file
testCallSitesMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryCallSitesMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline call_sites script detects a missing sample bam file
tryCallSitesMissingBamFileRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try cfsan_snp_pipeline call_sites with missing file
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1" &> "$logDir/callSites.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when the reads.sorted.deduped.indelrealigned.bam file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam does not exist"
    assertFileContains "$logDir/callSites.log" "Sample BAM file $tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam does not exist"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/callSites.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing sample bam file
testCallSitesMissingBamFileRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryCallSitesMissingBamFileRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing sample bam file
testCallSitesMissingBamFileRaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryCallSitesMissingBamFileRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing sample bam file
testCallSitesMissingBamFileRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryCallSitesMissingBamFileRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline map_reads script detects Samtools view failure.
tryMapReadsSamtoolsViewTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the SAM file and run cfsan_snp_pipeline map_reads
    echo "Garbage" > "$tempDir/samples/sample1/reads.sam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "samtools view"
    assertFileContains "$logDir/mapReads.log" "samtools view"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/mapReads.log" "Unsorted bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools view failure.
testMapReadsSamtoolsViewTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsSamtoolsViewTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools view failure.
testMapReadsSamtoolsViewTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsSamtoolsViewTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools view failure.
testMapReadsSamtoolsViewTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsSamtoolsViewTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects Samtools sort failure.
tryMapReadsSamtoolsSortTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the reads.unsorted.bam file and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.unsorted.bam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when reads.unsorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "samtools sort"
    assertFileContains "$logDir/mapReads.log" "samtools sort"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Sorted bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools sort failure.
testMapReadsSamtoolsSortTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsSamtoolsSortTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools sort failure.
testMapReadsSamtoolsSortTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsSamtoolsSortTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects Samtools sort failure.
testMapReadsSamtoolsSortTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsSamtoolsSortTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects Picard MarkDuplicates failure.
tryMapReadsPicardMarkDuplicatesTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the reads.sorted.bam file and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.bam"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when reads.sorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "picard"
    assertFileContains "$tempDir/error.log" "MarkDuplicates"
    assertFileContains "$logDir/mapReads.log" "picard"
    assertFileContains "$logDir/mapReads.log" "MarkDuplicates"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Deduped bam file is already freshly created "
}

# Verify the cfsan_snp_pipeline map_reads script detects Picard MarkDuplicates failure.
testMapReadsPicardMarkDuplicatesTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsPicardMarkDuplicatesTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects Picard MarkDuplicates failure.
testMapReadsPicardMarkDuplicatesTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsPicardMarkDuplicatesTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects Picard MarkDuplicates failure.
testMapReadsPicardMarkDuplicatesTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsPicardMarkDuplicatesTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run Picard MarkDuplicates.
tryMapReadsPicardMarkDuplicatesClasspathRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Clear CLASSPATH and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.bam"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when CLASSPATH unset." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute Picard. Define the path to Picard in the CLASSPATH environment variable."
    assertFileContains "$logDir/mapReads.log" "Error: cannot execute Picard. Define the path to Picard in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Deduped bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run PicardMarkDuplicates.
testMapReadsPicardMarkDuplicatesClasspathRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run PicardMarkDuplicates.
testMapReadsPicardMarkDuplicatesClasspathRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run PicardMarkDuplicates.
testMapReadsPicardMarkDuplicatesClasspathRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline map_reads script detects samtools index failure.
tryMapReadsSamtoolsIndexBamTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the reads.sorted.deduped.bam file and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when reads.sorted.deduped.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "samtools index"
    assertFileContains "$logDir/mapReads.log" "samtools index"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Bam file index is already freshly created "
}

# Verify the cfsan_snp_pipeline map_reads script detects samtools index failure.
testMapReadsSamtoolsIndexBamTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsSamtoolsIndexBamTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects samtools index failure.
testMapReadsSamtoolsIndexBamTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsSamtoolsIndexBamTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects samtools index failure.
testMapReadsSamtoolsIndexBamTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsSamtoolsIndexBamTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects RealignerTargetCreator failure.
tryMapReadsRealignerTargetCreatorTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the reads.sorted.deduped.bam file and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when reads.sorted.deduped.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "RealignerTargetCreator"
    assertFileContains "$logDir/mapReads.log" "RealignerTargetCreator"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Realign targets file is already freshly created "
}

# Verify the cfsan_snp_pipeline map_reads script detects RealignerTargetCreator failure.
testMapReadsRealignerTargetCreatorTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsRealignerTargetCreatorTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects RealignerTargetCreator failure.
testMapReadsRealignerTargetCreatorTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsRealignerTargetCreatorTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects RealignerTargetCreator failure.
testMapReadsRealignerTargetCreatorTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsRealignerTargetCreatorTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run RealignerTargetCreator.
tryMapReadsRealignerTargetCreatorClasspathRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Clear CLASSPATH and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when CLASSPATH unset." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute GATK RealignerTargetCreator. Define the path to GATK in the CLASSPATH environment variable."
    assertFileContains "$logDir/mapReads.log" "Error: cannot execute GATK RealignerTargetCreator. Define the path to GATK in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Realign targets file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run RealignerTargetCreator.
testMapReadsRealignerTargetCreatorClasspathRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsRealignerTargetCreatorClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run RealignerTargetCreator.
testMapReadsRealignerTargetCreatorClasspathRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsRealignerTargetCreatorClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run RealignerTargetCreator.
testMapReadsRealignerTargetCreatorClasspathRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsRealignerTargetCreatorClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline map_reads script detects IndelRealigner failure.
tryMapReadsIndelRealignerTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the realign.target.intervals file and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/realign.target.intervals"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when realign.target.intervals was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/mapReads.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "IndelRealigner"
    assertFileContains "$logDir/mapReads.log" "IndelRealigner"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Indelrealigned bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects IndelRealigner failure.
testMapReadsIndelRealignerTrapStop()
{
    export StopOnSampleError=true
    tryMapReadsIndelRealignerTrap 100
}

# Verify the cfsan_snp_pipeline map_reads script detects IndelRealigner failure.
testMapReadsIndelRealignerTrapNoStop()
{
    export StopOnSampleError=false
    tryMapReadsIndelRealignerTrap 98
}

# Verify the cfsan_snp_pipeline map_reads script detects IndelRealigner failure.
testMapReadsIndelRealignerTrapStopUnset()
{
    unset StopOnSampleError
    tryMapReadsIndelRealignerTrap 100
}


# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run IndelRealigner.
tryMapReadsIndelRealignerClasspathRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Clear CLASSPATH and run cfsan_snp_pipeline map_reads
    echo "Dummy" > "$tempDir/samples/sample1/reads.sam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/realign.target.intervals"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" &> "$logDir/mapReads.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline map_reads error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when CLASSPATH unset." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute GATK IndelRealigner. Define the path to GATK in the CLASSPATH environment variable."
    assertFileContains "$logDir/mapReads.log" "Error: cannot execute GATK IndelRealigner. Define the path to GATK in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/mapReads.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/mapReads.log" "Indelrealigned bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run IndelRealigner.
testMapReadsIndelRealignerClasspathRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMapReadsIndelRealignerClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run IndelRealigner.
testMapReadsIndelRealignerClasspathRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMapReadsIndelRealignerClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline map_reads script detects unset java classpath needed to run IndelRealigner.
testMapReadsIndelRealignerClasspathRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMapReadsIndelRealignerClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
tryCallSitesSamtoolsMpileupTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the reads.sorted.deduped.indelrealigned.bam file and run cfsan_snp_pipeline call_sites
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when reads.sorted.deduped.indelrealigned.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileNotContains "$logDir/callSites.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileContains "$tempDir/error.log" "samtools mpileup"
    assertFileContains "$logDir/callSites.log" "samtools mpileup"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/callSites.log" "Pileup file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testCallSitesSamtoolsMpileupTrapStop()
{
    export StopOnSampleError=true
    tryCallSitesSamtoolsMpileupTrap 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testCallSitesSamtoolsMpileupTrapNoStop()
{
    export StopOnSampleError=false
    tryCallSitesSamtoolsMpileupTrap 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testCallSitesSamtoolsMpileupTrapStopUnset()
{
    unset StopOnSampleError
    tryCallSitesSamtoolsMpileupTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
tryCallSitesVarscanRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Configure VarScan with invalid parameter settings and run cfsan_snp_pipeline call_sites
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    echo "Dummy" > "$tempDir/samples/sample1/reads.all.pileup"
    export VarscanMpileup2snp_ExtraParams="--min-coverage -1 --min-reads 99999999 --min-avg_qual -100 --min-var-freq 2 --output-vcf 2 --invalid-parameter"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"
    errorCode=$?
    export VarscanMpileup2snp_ExtraParams=""

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when varscan failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/callSites.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/callSites.log" "VarScan mpileup2snp"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/callSites.log" "Vcf file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testCallSitesVarscanRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryCallSitesVarscanRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testCallSitesVarscanRaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryCallSitesVarscanRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testCallSitesVarscanRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryCallSitesVarscanRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
tryCallSitesVarscanClasspathRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline call_sites with CLASSPATH unset
    echo "Dummy" > "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    echo "Dummy" > "$tempDir/samples/sample1/reads.all.pileup"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when CLASSPATH unset failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileContains "$logDir/callSites.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/callSites.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/callSites.log" "Vcf file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testCallSitesVarscanClasspathRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryCallSitesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testCallSitesVarscanClasspathRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryCallSitesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testCallSitesVarscanClasspathRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryCallSitesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
tryFilterRegionsPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Setup input files
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy vcf content" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample2/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample3/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample4/var.flt.vcf"

    # 1) var.flt_preserved.vcf =========
    # Make the output file var.flt_preserved.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_preserved.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_preserved.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_preserved.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_preserved.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for preserved SNPs"
    assertFileNotContains "$logDir/filterRegions.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_preserved.vcf"
    rm "$tempDir/error.log"

    # 2) var.flt_removed.vcf =========
    # Make the output file var.flt_removed.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_removed.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_removed.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_removed.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_removed.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for removed SNPs"
    assertFileNotContains "$logDir/filterRegions.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_removed.vcf"
    rm "$tempDir/error.log"
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testFilterRegionsPermissionTrapStop()
{
    export StopOnSampleError=true
    tryFilterRegionsPermissionTrap 100
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testFilterRegionsPermissionTrapNoStop()
{
    export StopOnSampleError=false

    expectErrorCode=0
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Setup input files
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy vcf content" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample2/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample3/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample4/var.flt.vcf"

    # 1) var.flt_preserved.vcf =========
    # Make the output file var.flt_preserved.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_preserved.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_preserved.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_preserved.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_preserved.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for preserved SNPs"
    assertFileNotContains "$logDir/filterRegions.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_preserved.vcf"
    rm "$tempDir/error.log"

    # 2) var.flt_removed.vcf =========
    # Make the output file var.flt_removed.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_removed.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_removed.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_removed.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_removed.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for removed SNPs"
    assertFileNotContains "$logDir/filterRegions.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_removed.vcf"
    rm "$tempDir/error.log"
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testFilterRegionsPermissionTrapStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsPermissionTrap 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
tryFilterRegionsMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline filter_regions with missing sampleDirectories.txt
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/filterRegions.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testFilterRegionsMissingSampleDirRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryFilterRegionsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testFilterRegionsMissingSampleDirRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryFilterRegionsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testFilterRegionsMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
tryFilterRegionsMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline filter_regions with missing reference
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy vcf content" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample2/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample3/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample4/var.flt.vcf"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/non-exist-reference" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/non-exist-reference does not exist"
    assertFileContains "$logDir/filterRegions.log" "Reference file $tempDir/non-exist-reference does not exist"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testFilterRegionsMissingReferenceRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryFilterRegionsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testFilterRegionsMissingReferenceRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryFilterRegionsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testFilterRegionsMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
tryFilterRegionsMissingOutgroupRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline filter_regions with missing outgroup samples file
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy vcf content" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample2/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample3/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample4/var.flt.vcf"
    cfsan_snp_pipeline filter_regions -g "$tempDir/outgroup" "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when file of outgroup samples file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "File of outgroup samples $tempDir/outgroup does not exist"
    assertFileContains "$logDir/filterRegions.log" "File of outgroup samples $tempDir/outgroup does not exist"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testFilterRegionsMissingOutgroupRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryFilterRegionsMissingOutgroupRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testFilterRegionsMissingOutgroupRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryFilterRegionsMissingOutgroupRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testFilterRegionsMissingOutgroupRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsMissingOutgroupRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
tryFilterRegionsMissingVcfRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline filter_regions -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when all var.flt.vcf were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$logDir/filterRegions.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testFilterRegionsMissingVcfRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryFilterRegionsMissingVcfRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testFilterRegionsMissingVcfRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryFilterRegionsMissingVcfRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testFilterRegionsMissingVcfRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsMissingVcfRaiseGlobalError 100
}


# Verify the filter_regions command detects missing some VCF files, but not all
tryFilterRegionsMissingVcfRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/mapReads.log"
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline filter_regions -- fail because of missing some, but not all VCF files. Only sample1 has been processed to this point.
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when some var.flt.vcf were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/filterRegions.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the filter_regions command detects missing some VCF files, but not all
testFilterRegionsMissingVcfRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryFilterRegionsMissingVcfRaiseSampleError 100
}

# Verify the filter_regions command detects missing some VCF files, but not all
testFilterRegionsMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline filter_regions -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/filterRegions.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when some var.flt.vcf were missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/filterRegions.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/filterRegions.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/filterRegions.log" "Use the -f option to force a rebuild"
}

# Verify the filter_regions command detects missing some VCF files, but not all
testFilterRegionsMissingVcfRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryFilterRegionsMissingVcfRaiseSampleError 100
}


# Verify the filter_regions command uses all the input vcf files to produce the outputs
# even when some of the samples are already fresh.
testFilterRegionsPartialRebuild()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> /dev/null

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Force timestamps to change so outputs are newer than inputs.
    # The files are small, quickly processed, and timestamps might not differ when we expect they will differ.
    touch -d '-12 day' $tempDir/reference/*.fasta
    touch -d '-11 day' $tempDir/reference/*.bt2
    touch -d '-10 day' $tempDir/samples/*/*.fastq
    touch -d  '-9 day' $tempDir/samples/*/reads.sam
    touch -d  '-8 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-7 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-6 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.deduped.bai
    touch -d  '-4 day' $tempDir/samples/*/reads.sorted.deduped.indelrealigned.bam
    touch -d  '-3 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-2 day' $tempDir/samples/*/var.flt.vcf
    touch -d  '-1 day' $tempDir/samples/*/var.flt_preserved.vcf
    touch -d  '-1 day' $tempDir/samples/*/var.flt_removed.vcf

    # Remove unwanted log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    rm -rf $logDir
    mkdir -p $logDir

    # Remove the results for one of the samples
    rm $tempDir/samples/sample1/var.flt_preserved.vcf

    # Re-run cfsan_snp_pipeline filter_regions -- this should only rebuild results for sample1, but it should use the var.flt.vcf input file for all samples
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" > "$logDir/filterRegions.log"

    # Verify log files
    verifyNonEmptyReadableFile "$logDir/filterRegions.log"
    assertFileNotContains "$logDir/filterRegions.log" "already freshly built"

    # Verify correct results
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/samples/sample1/var.flt_preserved.vcf" "$tempDir/expectedResults/samples/sample1/var.flt_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample2/var.flt_preserved.vcf" "$tempDir/expectedResults/samples/sample2/var.flt_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample3/var.flt_preserved.vcf" "$tempDir/expectedResults/samples/sample3/var.flt_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample4/var.flt_preserved.vcf" "$tempDir/expectedResults/samples/sample4/var.flt_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample1/var.flt_removed.vcf" "$tempDir/expectedResults/samples/sample1/var.flt_removed.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample2/var.flt_removed.vcf" "$tempDir/expectedResults/samples/sample2/var.flt_removed.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample3/var.flt_removed.vcf" "$tempDir/expectedResults/samples/sample3/var.flt_removed.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample4/var.flt_removed.vcf" "$tempDir/expectedResults/samples/sample4/var.flt_removed.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
}

# Verify the filter_regions command produces different results when one of the samples is an outgroup.
testFilterRegionsOutgroup()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> /dev/null

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Remember the output when there are no outgroups
    for d in $tempDir/samples/*; do
        verifyNonEmptyReadableFile $d/var.flt_preserved.vcf
        mv $d/var.flt_preserved.vcf $d/var.flt_preserved.vcf.save
    done

    # Force timestamps to change so outputs are newer than inputs.
    # The files are small, quickly processed, and timestamps might not differ when we expect they will differ.
    touch -d '-11 day' $tempDir/reference/*.fasta
    touch -d '-10 day' $tempDir/reference/*.bt2
    touch -d  '-9 day' $tempDir/samples/*/*.fastq
    touch -d  '-8 day' $tempDir/samples/*/reads.sam
    touch -d  '-7 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-6 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-4 day' $tempDir/samples/*/reads.sorted.deduped.bai
    touch -d  '-3 day' $tempDir/samples/*/reads.sorted.deduped.indelrealigned.bam
    touch -d  '-2 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-1 day' $tempDir/samples/*/var.flt.vcf

    # Remove unwanted log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    rm -rf $logDir
    mkdir -p $logDir

    # One of the samples is an outgroup
    outgroup="sample4" # this test only works when sample 4 is the outgroup
    echo $outgroup > "$tempDir/outgroup.txt"

    # Re-run cfsan_snp_pipeline filter_regions --
    cfsan_snp_pipeline filter_regions --out_group "$tempDir/outgroup.txt" "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" > "$logDir/filterRegions.log"

    # Verify log files
    verifyNonEmptyReadableFile "$logDir/filterRegions.log"
    assertFileNotContains "$logDir/filterRegions.log" "already freshly built"
    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"

    # Verify all preserved vcf files are different when there is an outgroup
    for d in $tempDir/samples/*; do
        verifyNonEmptyReadableFile $d/var.flt_preserved.vcf
        assertDifferentFiles $d/var.flt_preserved.vcf $d/var.flt_preserved.vcf.save --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    done

    # Verify the outgroup sample preserved vcf matches the original varscan output
    assertIdenticalFiles $tempDir/samples/$outgroup/var.flt_preserved.vcf $tempDir/samples/$outgroup/var.flt.vcf

    # Verify the outgroup sample removed vcf contains only a header
    grep -v "^#" $tempDir/samples/$outgroup/var.flt_removed.vcf > $tempDir/samples/$outgroup/var.flt_removed.noheader
    verifyEmptyFile $tempDir/samples/$outgroup/var.flt_removed.noheader
}


# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
tryMergeSitesPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Make the output file unwritable
    touch "$tempDir/snplist.txt"
    chmod -w "$tempDir/snplist.txt"

    # Try to create snplist.txt
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy vcf content" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample2/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample3/var.flt.vcf"
    echo "Dummy vcf content" > "$tempDir/samples/sample4/var.flt.vcf"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/mergeSites.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline merge_sites"
    assertFileNotContains "$logDir/mergeSites.log" "Error detected while running cfsan_snp_pipeline merge_sites"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/mergeSites.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snplist.txt"
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testMergeSitesPermissionTrapStop()
{
    export StopOnSampleError=true
    tryMergeSitesPermissionTrap 100
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testMergeSitesPermissionTrapNoStop()
{
    export StopOnSampleError=false
    tryMergeSitesPermissionTrap 100
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testMergeSitesPermissionTrapStopUnset()
{
    unset StopOnSampleError
    tryMergeSitesPermissionTrap 100
}



# Verify the merge_sites command detects failure.
tryMergeSitesMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_sites with missing sampleDirectories.txt
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/mergeSites.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/mergeSites.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/mergeSites.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingSampleDirRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMergeSitesMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingSampleDirRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMergeSitesMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMergeSitesMissingSampleDirRaiseGlobalError 100
}


# Verify the merge_sites command detects failure.
tryMergeSitesMissingVcfRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_sites -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/mergeSites.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$logDir/mergeSites.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/mergeSites.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMergeSitesMissingVcfRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMergeSitesMissingVcfRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMergeSitesMissingVcfRaiseGlobalError 100
}


# Verify the merge_sites command detects failure.
tryMergeSitesMissingVcfRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline merge_sites -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/mergeSites.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/mergeSites.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/mergeSites.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryMergeSitesMissingVcfRaiseSampleError 100
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline merge_sites -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/mergeSites.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/mergeSites.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/mergeSites.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/mergeSites.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testMergeSitesMissingVcfRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryMergeSitesMissingVcfRaiseSampleError 100
}


# Verify the call_consensus command traps on corrupt snplist file
tryCallConsensusCorruptSnplistTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline call_consensus -- fail because of corrupt snplist
    echo "Corrupt snplist content" > "$tempDir/snplist.txt"
    echo "Dummy pileup content" > "$tempDir/samples/sample1/reads.all.pileup"
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/callConsensus.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_consensus"
    assertFileNotContains "$logDir/callConsensus.log" "Error detected while running cfsan_snp_pipeline call_consensus"
    assertFileContains "$tempDir/error.log" "function read_snp_position_list at line"
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/callConsensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStop()
{
    export StopOnSampleError=true
    tryCallConsensusCorruptSnplistTrap 100
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapNoStop()
{
    export StopOnSampleError=false
    tryCallConsensusCorruptSnplistTrap 98
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStopUnset()
{
    unset StopOnSampleError
    tryCallConsensusCorruptSnplistTrap 100
}




# Verify the call_consensus command detects failure.
tryCallConsensusMissingSnpListRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline call_consensus -- fail because of missing snplist
    cfsan_snp_pipeline call_consensus -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/callConsensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the snplist file"
    assertFileContains "$logDir/callConsensus.log" "cannot call consensus without the snplist file"
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/callConsensus.log" "Snplist file snplist.txt does not exist"
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/callConsensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}


# Verify the call_consensus command does not fail merely because the snplist file is empty.
tryCallConsensusEmptySnpList()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline call_consensus with empty snplist
    touch "$tempDir/snplist.txt"
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/samples/sample1/consensus.fasta"  --vcfFileName "$tempDir/samples/sample1/consensus.vcf" "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/callConsensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when snplist was empty." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.vcf"
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/callConsensus.log" "cannot call consensus without the snplist file"
    assertFileNotContains "$logDir/callConsensus.log" "Snplist file $tempDir/snplist.txt is empty"
    assertFileContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/callConsensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStop()
{
    export StopOnSampleError=true
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListNoStop()
{
    export StopOnSampleError=false
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStopUnset()
{
    unset StopOnSampleError
    tryCallConsensusEmptySnpList 0
}


# Verify the call_consensus command detects missing pileup file.
tryCallConsensusMissingPileupRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    echo "fake snplist" > $tempDir/snplist.txt

    # Run cfsan_snp_pipeline call_consensus -- fail because of missing pileup
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/callConsensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when pileup file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the pileup file"
    assertFileContains "$logDir/callConsensus.log" "cannot call consensus without the pileup file"
    assertFileContains "$tempDir/error.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileContains "$logDir/callConsensus.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/callConsensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryCallConsensusMissingPileupRaiseSampleError 100
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryCallConsensusMissingPileupRaiseSampleError 98
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryCallConsensusMissingPileupRaiseSampleError 100
}


# Verify the call_consensus command detects missing exclude file.
tryCallConsensusMissingExcludeRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    echo "fake snplist" > $tempDir/snplist.txt
    mkdir -p "$tempDir/samples/sample1"
    echo "fake pileup" > "$tempDir/samples/sample1/reads.all.pileup"

    # Run cfsan_snp_pipeline call_consensus -- fail because of missing pileup
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline call_consensus -e "$tempDir/samples/sample1/excludeFile.vcf" -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/callConsensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when exclude file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the file of excluded positions"
    assertFileContains "$logDir/callConsensus.log" "cannot call consensus without the file of excluded positions"
    assertFileContains "$tempDir/error.log" "Exclude file $tempDir/samples/sample1/excludeFile.vcf does not exist"
    assertFileContains "$logDir/callConsensus.log" "Exclude file $tempDir/samples/sample1/excludeFile.vcf does not exist"
    assertFileNotContains "$logDir/callConsensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/callConsensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryCallConsensusMissingExcludeRaiseSampleError 100
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryCallConsensusMissingExcludeRaiseSampleError 98
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryCallConsensusMissingExcludeRaiseSampleError 100
}


# Verify the merge_vcfs command detects failure.
tryMergeVcfsCorruptVcfTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with corrupt consensus.vcf
    sleep 1
    sed -i 's/1/@@@/g' "$tempDir/samples/sample1/consensus.vcf"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when consensus.vcf was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline merge_vcfs."
    assertFileNotContains "$logDir/mergeVcfs.log" "Error detected while running cfsan_snp_pipeline merge_vcfs."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsCorruptVcfTrapStop()
{
    export StopOnSampleError=true
    tryMergeVcfsCorruptVcfTrap 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsCorruptVcfTrapNoStop()
{
    export StopOnSampleError=false
    tryMergeVcfsCorruptVcfTrap 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsCorruptVcfTrapStopUnset()
{
    unset StopOnSampleError
    tryMergeVcfsCorruptVcfTrap 100
}


# Verify the merge_vcfs command detects failure.
tryMergeVcfsMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with missing sampleDirectories.txt
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist."
    assertFileContains "$logDir/mergeVcfs.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsMissingSampleDirRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryMergeVcfsMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsMissingSampleDirRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMergeVcfsMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfsMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryMergeVcfsMissingSampleDirRaiseGlobalError 100
}


# Verify the merge_vcfs command detects a missing consensus VCF file.
tryMergeVcfsMissingVcfRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with missing vcf files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when vcf file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcfs.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects a missing consensus VCF file.
testMergeVcfsMissingVcfRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryMergeVcfsMissingVcfRaiseSampleError 100
}

# Verify the merge_vcfs command detects a missing consensus VCF file - but continues running.
testMergeVcfsMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with missing vcf files
    rm "$tempDir/samples/sample1/consensus.vcf"
    rm "$tempDir/snpma.vcf"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command keeps running when only one vcf file is missing
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when vcf file was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcfs.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
}

# Verify the merge_vcfs command detects a missing consensus VCF file.
testMergeVcfsMissingVcfRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryMergeVcfsMissingVcfRaiseSampleError 100
}


# Verify the merge_vcfs command simply copies the input consensus VCF file when there is only one sample
testMergeVcfsOnlyOneSample()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export StopOnSampleError=true
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with only one vcf file
    mkdir -p "$tempDir/samples/sample1"
    echo "$tempDir/samples/sample1" > "$tempDir/sampleDirectories.txt"
    echo "Dummy VCF contents" > "$tempDir/samples/sample1/consensus.vcf"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command copies the input consensus VCF file when there is only one sample
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when there was only one vcf file." 0 $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcfs.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    assertIdenticalFiles "$tempDir/samples/sample1/consensus.vcf" "$tempDir/snpma.vcf"
}


# Verify the merge_vcfs command detects all the consensus VCF files missing
tryMergeVcfsZeroGoodSamplesRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with no consensus vcf files
    echo   "$tempDir/samples/sample1" > "$tempDir/sampleDirectories.txt"
    echo   "$tempDir/samples/sample2" >> "$tempDir/sampleDirectories.txt"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcfs.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when no good VCF files." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "There are no vcf files to merge."
    assertFileContains "$logDir/mergeVcfs.log" "There are no vcf files to merge."
    assertFileNotContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcfs.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
#testMergeVcfsZeroGoodSamplesRaiseGlobalErrorStop()
#{
#    # Nothing to test, StopOnSampleError must be false to test this code path.
#    # Otherwise, the first missing VCF file will trigger stop upon sample error -- already tested.
#}

# Verify the merge_vcfs command detects failure.
testMergeVcfsZeroGoodSamplesRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryMergeVcfsZeroGoodSamplesRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
#testMergeVcfsZeroGoodSamplesRaiseGlobalErrorStopUnset()
#{
#    # Nothing to test, StopOnSampleError must be false to test this code path.
#    # Otherwise, the first missing VCF file will trigger stop upon sample error -- already tested.
#}



# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
trySnpMatrixPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Make the output file unwritable
    touch "$tempDir/snpma.fasta"
    chmod -w "$tempDir/snpma.fasta"

    # Try to create snplist.txt
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy content" > "$tempDir/samples/sample1/consensus.fasta"
    echo "Dummy content" > "$tempDir/samples/sample2/consensus.fasta"
    echo "Dummy content" > "$tempDir/samples/sample3/consensus.fasta"
    echo "Dummy content" > "$tempDir/samples/sample4/consensus.fasta"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/snpMatrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline snp_matrix"
    assertFileNotContains "$logDir/snpMatrix.log" "Error detected while running cfsan_snp_pipeline snp_matrix"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/snpMatrix.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snpma.fasta"
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testSnpMatrixPermissionTrapStop()
{
    export StopOnSampleError=true
    trySnpMatrixPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testSnpMatrixPermissionTrapNoStop()
{
    export StopOnSampleError=false
    trySnpMatrixPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testSnpMatrixPermissionTrapStopUnset()
{
    unset StopOnSampleError
    trySnpMatrixPermissionTrap 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
trySnpMatrixMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline snp_matrix with missing sampleDirectories.txt
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/snpMatrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/snpMatrix.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/snpMatrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingSampleDirRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    trySnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingSampleDirRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    trySnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    trySnpMatrixMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
trySnpMatrixMissingConsensusRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    #cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/indexRef.log"
    #cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    #cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/callSites.log"

    # Run cfsan_snp_pipeline snp_matrix -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpmap.fasta"  "$tempDir/sampleDirList.txt" &> "$logDir/snpMatrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when all consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample2/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample3/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample2/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample3/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: all 4 consensus fasta files were missing or empty"
    assertFileContains "$logDir/snpMatrix.log" "Error: all 4 consensus fasta files were missing or empty"
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/snpMatrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    trySnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    trySnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    trySnpMatrixMissingConsensusRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
trySnpMatrixMissingConsensusRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Run cfsan_snp_pipeline snp_matrix -- fail because of missing some, but not all consensus fasta files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    rm "$tempDir/samples/sample1/consensus.fasta"
    rm "$tempDir/samples/sample4/consensus.fasta"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/snpMatrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when some consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/snpMatrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/snpMatrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    trySnpMatrixMissingConsensusRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Run cfsan_snp_pipeline snp_matrix -- fail because of missing some, but not all consensus fasta files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    rm "$tempDir/samples/sample1/consensus.fasta"
    rm "$tempDir/samples/sample4/consensus.fasta"
    rm "$tempDir/snpma.fasta"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/snpMatrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when some consensus.fasta missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/snpMatrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/snpMatrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/snpMatrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testSnpMatrixMissingConsensusRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    trySnpMatrixMissingConsensusRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
trySnpReferencePermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Make the output file unwritable
    touch "$tempDir/referenceSNP.fasta"
    chmod -w "$tempDir/referenceSNP.fasta"

    # Try to create snplist.txt
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy content" > "$tempDir/snplist.txt"
    cfsan_snp_pipeline snp_reference -l "$tempDir/snplist.txt" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snpReference.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline snp_reference"
    assertFileNotContains "$logDir/snpReference.log" "Error detected while running cfsan_snp_pipeline snp_reference"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/snpReference.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/referenceSNP.fasta"
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testSnpReferencePermissionTrapStop()
{
    export StopOnSampleError=true
    trySnpReferencePermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testSnpReferencePermissionTrapNoStop()
{
    export StopOnSampleError=false
    trySnpReferencePermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testSnpReferencePermissionTrapStopUnset()
{
    unset StopOnSampleError
    trySnpReferencePermissionTrap 100
}


# Verify the cfsan_snp_pipeline snp_reference script detects failure.
trySnpReferenceMissingSnpListRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline snp_reference -- fail because of missing snplist
    cfsan_snp_pipeline snp_reference -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snpReference.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileNotContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/snpReference.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the snplist file"
    assertFileContains "$logDir/snpReference.log" "cannot create the snp reference sequence without the snplist file"
    assertFileNotContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/snpReference.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingSnpListRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    trySnpReferenceMissingSnpListRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingSnpListRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    trySnpReferenceMissingSnpListRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    trySnpReferenceMissingSnpListRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_reference script detects failure.
trySnpReferenceMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline snp_reference -- fail because of missing reference
    echo "Dummy snplist content" > "$tempDir/snplist"
    rm "$tempDir/reference/lambda_virus.fasta"
    cfsan_snp_pipeline snp_reference -l "$tempDir/snplist" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snpReference.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileNotContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/snpReference.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileContains "$logDir/snpReference.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileNotContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/snpReference.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingReferenceRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    trySnpReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingReferenceRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    trySnpReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testSnpReferenceMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    trySnpReferenceMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
tryCollectMetricsMissingSampleDirRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately remove directory
    rm -rf "$tempDir/samples/sample1"

    # Try to collect metrics
    echo "Dummy consensus.fasta content" > "$tempDir/consensus.fasta"
    cfsan_snp_pipeline collect_metrics -c "$tempDir/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metricsn error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the sample directory was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileNotContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileContains "$logDir/collectMetrics.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileNotContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectMetricsMissingSampleDirRaiseSampleErrorStop()
{
    export StopOnSampleError=true
    tryCollectMetricsMissingSampleDirRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectMetricsMissingSampleDirRaiseSampleErrorNoStop()
{
    export StopOnSampleError=false
    tryCollectMetricsMissingSampleDirRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectMetricsMissingSampleDirRaiseSampleErrorStopUnset()
{
    unset StopOnSampleError
    tryCollectMetricsMissingSampleDirRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
tryCollectMetricsMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately remove reference
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to collect metrics
    echo "Dummy consensus.fasta content" > "$tempDir/samples/sample1/consensus.fasta"
    cfsan_snp_pipeline collect_metrics -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileNotContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/collectMetrics.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectMetricsMissingReferenceRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryCollectMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectMetricsMissingReferenceRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryCollectMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectMetricsMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryCollectMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
tryCollectMetricsMissingInputFiles()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Try to collect metrics
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."

    assertFileContains "$logDir/collectMetrics.log" "SAM file reads.sam was not found"
    assertFileContains "$logDir/collectMetrics.log" "Deduped BAM file reads.sorted.deduped.bam was not found"
    assertFileContains "$logDir/collectMetrics.log" "BAM file reads.sorted.bam was not found"
    assertFileContains "$logDir/collectMetrics.log" "Pileup file reads.all.pileup was not found"
    assertFileContains "$logDir/collectMetrics.log" "VCF file var.flt.vcf was not found"
    assertFileContains "$logDir/collectMetrics.log" "VCF file var.flt_preserved.vcf was not found"
    assertFileContains "$logDir/collectMetrics.log" "Consensus VCF file consensus.vcf was not found"
    assertFileContains "$logDir/collectMetrics.log" "Consensus VCF file consensus_preserved.vcf was not found"
    assertFileContains "$logDir/collectMetrics.log" "Consensus fasta file consensus.fasta was not found"
    assertFileContains "$logDir/collectMetrics.log" "Consensus fasta file consensus_preserved.fasta was not found"

    assertFileContains "$tempDir/samples/sample1/metrics" "SAM file reads.sam was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Deduped BAM file reads.sorted.deduped.bam was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "BAM file reads.sorted.bam was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Pileup file reads.all.pileup was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "VCF file var.flt.vcf was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "VCF file var.flt_preserved.vcf was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus VCF file consensus.vcf was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus VCF file consensus_preserved.vcf was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus fasta file consensus.fasta was not found"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus fasta file consensus_preserved.fasta was not found"

    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSamplePreserved=$"

    assertFileNotContains "$tempDir/samples/sample1/metrics" "Cannot calculate"
    assertFileNotContains "$logDir/collectMetrics.log" "Cannot calculate"
    assertFileContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectMetricsMissingInputFilesStop()
{
    export StopOnSampleError=true
    tryCollectMetricsMissingInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectMetricsMissingInputFilesNoStop()
{
    export StopOnSampleError=false
    tryCollectMetricsMissingInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectMetricsMissingInputFilesStopUnset()
{
    unset StopOnSampleError
    tryCollectMetricsMissingInputFiles 0
}


# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
tryCollectMetricsEmptyInputFiles()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately create empty files
    touch "$tempDir/samples/sample1/reads.sam"
    touch "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    touch "$tempDir/samples/sample1/reads.all.pileup"
    touch "$tempDir/samples/sample1/var.flt.vcf"
    touch "$tempDir/samples/sample1/var.flt_preserved.vcf"
    touch "$tempDir/samples/sample1/consensus.vcf"
    touch "$tempDir/samples/sample1/consensus_preserved.vcf"
    touch "$tempDir/samples/sample1/consensus.fasta"
    touch "$tempDir/samples/sample1/consensus_preserved.fasta"

    # Try to collect metrics
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."

    assertFileContains "$logDir/collectMetrics.log" "SAM file reads.sam is empty"
    assertFileContains "$logDir/collectMetrics.log" "Deduped BAM file reads.sorted.deduped.bam is empty"
    assertFileContains "$logDir/collectMetrics.log" "BAM file reads.sorted.bam is empty"
    assertFileContains "$logDir/collectMetrics.log" "Pileup file reads.all.pileup is empty"
    assertFileContains "$logDir/collectMetrics.log" "VCF file var.flt.vcf is empty"
    assertFileContains "$logDir/collectMetrics.log" "VCF file var.flt_preserved.vcf is empty"
    assertFileContains "$logDir/collectMetrics.log" "Consensus VCF file consensus.vcf is empty"
    assertFileContains "$logDir/collectMetrics.log" "Consensus VCF file consensus_preserved.vcf is empty"
    assertFileContains "$logDir/collectMetrics.log" "Consensus fasta file consensus.fasta is empty"
    assertFileContains "$logDir/collectMetrics.log" "Consensus fasta file consensus_preserved.fasta is empty"

    assertFileContains "$tempDir/samples/sample1/metrics" "SAM file reads.sam is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Deduped BAM file reads.sorted.deduped.bam is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "BAM file reads.sorted.bam is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Pileup file reads.all.pileup is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "VCF file var.flt.vcf is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "VCF file var.flt_preserved.vcf is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus VCF file consensus.vcf is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus VCF file consensus_preserved.vcf is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus fasta file consensus.fasta is empty"
    assertFileContains "$tempDir/samples/sample1/metrics" "Consensus fasta file consensus_preserved.fasta is empty"

    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSamplePreserved=$"

    assertFileNotContains "$tempDir/samples/sample1/metrics" "Cannot calculate"
    assertFileNotContains "$logDir/collectMetrics.log" "Cannot calculate"
    assertFileContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectMetricsEmptyInputFilesStop()
{
    export StopOnSampleError=true
    tryCollectMetricsEmptyInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectMetricsEmptyInputFilesNoStop()
{
    export StopOnSampleError=false
    tryCollectMetricsEmptyInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectMetricsEmptyInputFilesStopUnset()
{
    unset StopOnSampleError
    tryCollectMetricsEmptyInputFiles 0
}


# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
tryCollectMetricsCorruptInputFiles()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt files
    echo "Garbage" > "$tempDir/samples/sample1/reads.sam"
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.bam"
    echo "Garbage" > "$tempDir/samples/sample1/reads.all.pileup"
    echo "Garbage" > "$tempDir/samples/sample1/var.flt.vcf"
    echo "Garbage" > "$tempDir/samples/sample1/var.flt_preserved.vcf"
    echo "Garbage" > "$tempDir/samples/sample1/consensus.vcf"
    echo "Garbage" > "$tempDir/samples/sample1/consensus_preserved.vcf"
    echo "Garbage" > "$tempDir/samples/sample1/consensus.fasta"
    echo "Garbage" > "$tempDir/samples/sample1/consensus_preserved.fasta"

    # Try to collect metrics
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."
    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=0$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSamplePreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "Cannot calculate number of reads and %mapped"
    assertFileContains "$logDir/collectMetrics.log" "Cannot calculate number of reads and %mapped"
    assertFileContains "$tempDir/samples/sample1/metrics" "Cannot calculate mean insert size"
    assertFileContains "$logDir/collectMetrics.log" "Cannot calculate mean insert size"
    assertFileContains "$tempDir/samples/sample1/metrics" "Cannot calculate mean pileup depth"
    assertFileContains "$logDir/collectMetrics.log" "Cannot calculate mean pileup depth"
    assertFileContains "$logDir/collectMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectMetricsCorruptInputFilesStop()
{
    export StopOnSampleError=true
    tryCollectMetricsCorruptInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectMetricsCorruptInputFilesNoStop()
{
    export StopOnSampleError=false
    tryCollectMetricsCorruptInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectMetricsCorruptInputFilesStopUnset()
{
    unset StopOnSampleError
    tryCollectMetricsCorruptInputFiles 0
}


# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
tryCombineMetricsMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline combine_metrics with missing sampleDirectories.txt
    cfsan_snp_pipeline combine_metrics "$tempDir/sampleDirectories.txt" &> "$logDir/combineMetrics.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics failed."
    assertFileNotContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/combineMetrics.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineMetricsMissingSampleDirRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryCombineMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineMetricsMissingSampleDirRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryCombineMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineMetricsMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryCombineMetricsMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
tryCombineMetricsMissingSampleMetricsRaiseSampleWarning()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Try to combine metrics
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    touch $tempDir/samples/sample4/metrics
    cfsan_snp_pipeline combine_metrics -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline combine_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when the sample metrics file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics warning"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics failed"
    assertFileNotContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics failed"
    assertFileContains "$tempDir/metrics.tsv" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$tempDir/error.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$logDir/combineMetrics.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$tempDir/metrics.tsv" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$tempDir/error.log" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$logDir/combineMetrics.log" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineMetricsMissingSampleMetricsRaiseSampleWarningStop()
{
    export StopOnSampleError=true
    tryCombineMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineMetricsMissingSampleMetricsRaiseSampleWarningNoStop()
{
    export StopOnSampleError=false
    tryCombineMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineMetricsMissingSampleMetricsRaiseSampleWarningStopUnset()
{
    unset StopOnSampleError
    tryCombineMetricsMissingSampleMetricsRaiseSampleWarning 0
}


# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
tryCombineMetricsPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Make the output file unwritable
    touch "$tempDir/metrics.tsv"
    chmod -w "$tempDir/metrics.tsv"

    # Try to combine metrics
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy snpma.fasta content" > "$tempDir/snpma.fasta"
    cfsan_snp_pipeline combine_metrics -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline combine_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline combine_metrics"
    assertFileNotContains "$logDir/combineMetrics.log" "Error detected while running cfsan_snp_pipeline combine_metrics"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineMetrics.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/metrics.tsv"
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineMetricsPermissionTrapStop()
{
    export StopOnSampleError=true
    tryCombineMetricsPermissionTrap 100
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineMetricsPermissionTrapNoStop()
{
    export StopOnSampleError=false
    tryCombineMetricsPermissionTrap 100
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineMetricsPermissionTrapStopUnset()
{
    unset StopOnSampleError
    tryCombineMetricsPermissionTrap 100
}


# Verify the cfsan_snp_pipeline distance script detects missing input file
tryDistanceMissingInputRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline distance with missing snpma.fasta
    cfsan_snp_pipeline distance -p pp -m mm "$tempDir/snpma.fasta" &> "$logDir/distance.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when input snp matrix file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline distance failed."
    assertFileNotContains "$logDir/distance.log" "cfsan_snp_pipeline distance failed"
    assertFileContains "$tempDir/error.log" "Error: cannot calculate sequence distances without the snp matrix file"
    assertFileContains "$logDir/distance.log" "Error: cannot calculate sequence distances without the snp matrix file"
    assertFileContains "$tempDir/error.log" "SNP matrix file $tempDir/snpma.fasta does not exist"
    assertFileContains "$logDir/distance.log" "SNP matrix file $tempDir/snpma.fasta does not exist"
    assertFileNotContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/distance.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testDistanceMissingInputRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryDistanceMissingInputRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testDistanceMissingInputRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryDistanceMissingInputRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testDistanceMissingInputRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryDistanceMissingInputRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline distance script detects no output file options
tryDistanceMissingOutputOptionsRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Create dummpy input file
    touch "$tempDir/snpma.fasta"

    # Run cfsan_snp_pipeline distance with no output options
    cfsan_snp_pipeline distance "$tempDir/snpma.fasta" &> "$logDir/distance.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when both output options were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline distance failed."
    assertFileNotContains "$logDir/distance.log" "cfsan_snp_pipeline distance failed"
    assertFileContains "$tempDir/error.log" "Error: no output file specified"
    assertFileContains "$logDir/distance.log" "Error: no output file specified"
    assertFileNotContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/distance.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testDistanceMissingOutputOptionsRaiseGlobalErrorStop()
{
    export StopOnSampleError=true
    tryDistanceMissingOutputOptionsRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testDistanceMissingOutputOptionsRaiseGlobalErrorNoStop()
{
    export StopOnSampleError=false
    tryDistanceMissingOutputOptionsRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testDistanceMissingOutputOptionsRaiseGlobalErrorStopUnset()
{
    unset StopOnSampleError
    tryDistanceMissingOutputOptionsRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
tryDistancePermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Make the output file unwritable
    touch "$tempDir/pairwise"
    chmod -w "$tempDir/pairwise"

    # Try to create snp distances
    echo "> Sequence" > "$tempDir/snpma.fasta"
    echo "ACGT" >> "$tempDir/snpma.fasta"
    cfsan_snp_pipeline distance -p "$tempDir/pairwise" "$tempDir/snpma.fasta" &> "$logDir/distance.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline distance error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline distance"
    assertFileNotContains "$logDir/distance.log" "Error detected while running cfsan_snp_pipeline distance"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileContains "$logDir/distance.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/distance.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/pairwise"
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testDistancePermissionTrapStop()
{
    export StopOnSampleError=true
    tryDistancePermissionTrap 100
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testDistancePermissionTrapNoStop()
{
    export StopOnSampleError=false
    tryDistancePermissionTrap 100
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testDistancePermissionTrapStopUnset()
{
    unset StopOnSampleError
    tryDistancePermissionTrap 100
}


# Verify the run_snp_pipeline.sh script detects missing reference
tryRunSnpPipelineMissingReferenceRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Deliberately try run with missing file
    rm "$tempDir/reference/lambda_virus.fasta"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing configuration file
tryRunSnpPipelineMissingConfigurationFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/not-exist.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the configuration file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "configuration file $tempDir/not-exist.conf does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "configuration file $tempDir/not-exist.conf does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir
    cfsan_snp_pipeline data configurationFile $tempDir

    # Misconfigure the aligner
    echo "SnpPipeline_Aligner=garbage" >> "$tempDir/snppipeline.conf"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the aligner was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Config file error in SnpPipeline_Aligner parameter"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Config file error in SnpPipeline_Aligner parameter"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing file of sample directories
tryRunSnpPipelineMissingSampleDirFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -S "not-exist-file" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the file of sample directories was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "file of samples directories, not-exist-file, does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "file of samples directories, not-exist-file, does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script validates the file of sample directories
tryRunSnpPipelineValidateSampleDirFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Create some error conditions
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    rm -rf "$tempDir/samples/sample1"  # sample1 directory does not exist
    rm $tempDir/samples/sample2/*.*  # sample2 directory is empty
    rm $tempDir/samples/sample3/*.*  # sample3 directory does not contain any fastq files
    echo "dummy" > "$tempDir/samples/sample3/not-a-fastq-file"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -S "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the file of sample directories failed validation." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample1 does not exist."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample1 does not exist."
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample2 is empty."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample2 is empty."
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample3 does not contain any fastq files."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample3 does not contain any fastq files."
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineValidateSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir
    cp -r "$tempDir/samples/sample1" "$tempDir/samples/sample5"  # extra sample so at least 2 will be merged

    # Create some error conditions
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    rm -rf "$tempDir/samples/sample1"  # sample1 directory does not exist
    rm $tempDir/samples/sample2/*.*  # sample2 directory is empty
    rm $tempDir/samples/sample3/*.*  # sample3 directory does not contain any fastq files
    echo "dummy" > "$tempDir/samples/sample3/not-a-fastq-file"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -S "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the file of sample directories failed validation." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample1 does not exist."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample1 does not exist."
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample2 is empty."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample2 is empty."
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample3 does not contain any fastq files."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Sample directory $tempDir/samples/sample3 does not contain any fastq files."
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample2/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample3/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample5/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/consensus.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample5/consensus.vcf"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    verifyNonEmptyReadableFile "$tempDir/metrics.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix.tsv"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/consensus_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample5/consensus_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/consensus_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample5/consensus_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise_preserved.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix_preserved.tsv"

    assertFileContains "$tempDir/snpma.fasta" "sample4"
    assertFileContains "$tempDir/snpma.fasta" "sample5"
    assertFileContains "$tempDir/snpma.vcf" "sample4"
    assertFileContains "$tempDir/snpma.vcf" "sample5"

    assertFileContains "$tempDir/snpma_preserved.fasta" "sample4"
    assertFileContains "$tempDir/snpma_preserved.fasta" "sample5"
    assertFileContains "$tempDir/snpma_preserved.vcf" "sample4"
    assertFileContains "$tempDir/snpma_preserved.vcf" "sample5"

    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "See the log file $tempDir/error.log for a summary of errors."
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineValidateSampleDirFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing directory of samples
tryRunSnpPipelineMissingSamplesDirRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "not-exist-dir" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the directory of samples was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Samples directory not-exist-dir does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Samples directory not-exist-dir does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects empty directory of samples
tryRunSnpPipelineEmptySamplesDirRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    mkdir -p "$tempDir/emptySamplesDir"
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/emptySamplesDir" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the directory of samples was empty." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Samples directory $tempDir/emptySamplesDir is empty"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Samples directory $tempDir/emptySamplesDir is empty"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    find $tempDir -name *.fastq -delete
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh returned incorrect error code when the directory of samples was empty." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/error.log" "run_snp_pipeline.sh failed"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "run_snp_pipeline.sh failed"
    assertFileContains "$tempDir/error.log" "Samples directory $tempDir/samples does not contain subdirectories with fastq files"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Samples directory $tempDir/samples does not contain subdirectories with fastq files"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}


# Verify run_snp_pipeline.sh trap handling
tryRunSnpPipelineTrapPrepReferenceTrap()
{
    expectErrorCode=$1

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when the input fasta was corrupt." $expectErrorCode $errorCode

    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/error.log" "error|failed"
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "error|failed"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2-build"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"

    assertFileContains "$tempDir/error.log" "See also the log files in directory $tempDir/logs"
    assertFileContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/reference/lambda_virus.1.bt2"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample2/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample3/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample4/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"

    verifyNonExistingFile "$tempDir/samples/sample1/var.flt_preserved.vcf"
    verifyNonExistingFile "$tempDir/snplist_preserved.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.vcf"
    verifyNonExistingFile "$tempDir/snpma_preserved.fasta"
    verifyNonExistingFile "$tempDir/snpma_preserved.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP_preserved.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}


# Verify run_snp_pipeline.sh trap handling
tryRunSnpPipelineTrapAlignSampleToReferenceTrap()
{
    expectErrorCode=$1

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Deliberately corrupt the fastq files
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample2/sample2_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample4/sample4_1.fastq"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when the input fastq files were corrupt." $expectErrorCode $errorCode
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline call_sites"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline merge_sites"

    assertFileContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"

    verifyNonExistingFile "$tempDir/samples/sample1/var.flt_preserved.vcf"
    verifyNonExistingFile "$tempDir/snplist_preserved.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.vcf"
    verifyNonExistingFile "$tempDir/snpma_preserved.fasta"
    verifyNonExistingFile "$tempDir/snpma_preserved.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP_preserved.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapAlignSampleToReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapNoStopAllFail()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    expectErrorCode=1

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Deliberately corrupt the fastq files
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample2/sample2_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample4/sample4_1.fastq"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when all the input fastq files were corrupt." $expectErrorCode $errorCode
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline call_sites"

    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample2/reads.sorted.deduped.indelrealigned.bam"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample3/reads.sorted.deduped.indelrealigned.bam"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample4/reads.sorted.deduped.indelrealigned.bam"

    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed|cfsan_snp_pipeline merge_sites failed"  # either/or
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: all 4 VCF files were missing or empty"

    assertFileContains "$tempDir/error.log" "See also the log files in directory $tempDir/logs"
    assertFileContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bai"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"

    verifyNonExistingFile "$tempDir/samples/sample1/var.flt_preserved.vcf"
    verifyNonExistingFile "$tempDir/snplist_preserved.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.vcf"
    verifyNonExistingFile "$tempDir/snpma_preserved.fasta"
    verifyNonExistingFile "$tempDir/snpma_preserved.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP_preserved.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapNoStopSomeFail()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    expectErrorCode=0

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Deliberately corrupt some of the fastq files
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample4/sample4_1.fastq"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when some the input fastq files were corrupt." $expectErrorCode $errorCode
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline call_sites"

    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample1/reads.sorted.deduped.indelrealigned.bam"
    assertFileNotContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample2/reads.sorted.deduped.indelrealigned.bam"
    assertFileNotContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample3/reads.sorted.deduped.indelrealigned.bam"
    assertFileContains "$tempDir/error.log" "Sample BAM file $tempDir/samples/sample4/reads.sorted.deduped.indelrealigned.bam"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline call_sites finished"

    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 VCF files were missing or empty"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline filter_regions finished"

    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 VCF files were missing or empty"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline merge_sites finished"

    assertFileNotContains "$tempDir/error.log" "See also the log files in directory $tempDir/logs"
    assertFileNotContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample4/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample4/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/consensus.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/consensus.vcf"
    verifyNonExistingFile "$tempDir/samples/sample4/consensus.vcf"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    verifyNonEmptyReadableFile "$tempDir/metrics.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix.tsv"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/consensus_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/consensus_preserved.fasta"
    verifyNonExistingFile "$tempDir/samples/sample4/consensus_preserved.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/consensus_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/consensus_preserved.vcf"
    verifyNonExistingFile "$tempDir/samples/sample4/consensus_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_pairwise_preserved.tsv"
    verifyNonEmptyReadableFile "$tempDir/snp_distance_matrix_preserved.tsv"

    assertFileContains "$tempDir/snpma.fasta" "sample2"
    assertFileContains "$tempDir/snpma.fasta" "sample3"
    assertFileContains "$tempDir/snpma.vcf" "sample2"
    assertFileContains "$tempDir/snpma.vcf" "sample3"

    assertFileContains "$tempDir/snpma_preserved.fasta" "sample2"
    assertFileContains "$tempDir/snpma_preserved.fasta" "sample3"
    assertFileContains "$tempDir/snpma_preserved.vcf" "sample2"
    assertFileContains "$tempDir/snpma_preserved.vcf" "sample3"

    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "See the log file $tempDir/error.log for a summary of errors."
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapAlignSampleToReferenceTrap 1
}


# Verify run_snp_pipeline generates correct results for the lambda data set
testRunSnpPipelineLambda()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Simulate left-over errors from a prior run
    echo "Old errors from prior run" > "$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify correct mirroring
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample1/sample1_2.fastq" "$tempDir/samples/sample1/sample1_2.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample2/sample2_1.fastq" "$tempDir/samples/sample2/sample2_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample2/sample2_2.fastq" "$tempDir/samples/sample2/sample2_2.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample3/sample3_1.fastq" "$tempDir/samples/sample3/sample3_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample3/sample3_2.fastq" "$tempDir/samples/sample3/sample3_2.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample4/sample4_1.fastq" "$tempDir/samples/sample4/sample4_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample4/sample4_2.fastq" "$tempDir/samples/sample4/sample4_2.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/reference/lambda_virus.fasta"    "$tempDir/reference/lambda_virus.fasta"

    # Verify correct results
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snplist.txt"                   "$tempDir/expectedResults/snplist.txt"
    assertIdenticalFiles "$tempDir/snpma.fasta"                   "$tempDir/expectedResults/snpma.fasta"
    assertIdenticalFiles "$tempDir/snpma.vcf"                     "$tempDir/expectedResults/snpma.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source --ignore-matching-lines=##bcftools
    assertIdenticalFiles "$tempDir/referenceSNP.fasta"            "$tempDir/expectedResults/referenceSNP.fasta"
    assertIdenticalFiles "$tempDir/metrics.tsv"                   "$tempDir/expectedResults/metrics.tsv"
    assertIdenticalFiles "$tempDir/snp_distance_pairwise.tsv"     "$tempDir/expectedResults/snp_distance_pairwise.tsv"
    assertIdenticalFiles "$tempDir/snp_distance_matrix.tsv"       "$tempDir/expectedResults/snp_distance_matrix.tsv"

    assertIdenticalFiles "$tempDir/samples/sample1/reads.sam"     "$tempDir/expectedResults/samples/sample1/reads.sam" --ignore-matching-lines=bowtie
    assertIdenticalFiles "$tempDir/samples/sample2/reads.sam"     "$tempDir/expectedResults/samples/sample2/reads.sam" --ignore-matching-lines=bowtie
    assertIdenticalFiles "$tempDir/samples/sample3/reads.sam"     "$tempDir/expectedResults/samples/sample3/reads.sam" --ignore-matching-lines=bowtie
    assertIdenticalFiles "$tempDir/samples/sample4/reads.sam"     "$tempDir/expectedResults/samples/sample4/reads.sam" --ignore-matching-lines=bowtie

    assertIdenticalFiles "$tempDir/samples/sample1/var.flt.vcf"   "$tempDir/expectedResults/samples/sample1/var.flt.vcf"
    assertIdenticalFiles "$tempDir/samples/sample2/var.flt.vcf"   "$tempDir/expectedResults/samples/sample2/var.flt.vcf"
    assertIdenticalFiles "$tempDir/samples/sample3/var.flt.vcf"   "$tempDir/expectedResults/samples/sample3/var.flt.vcf"
    assertIdenticalFiles "$tempDir/samples/sample4/var.flt.vcf"   "$tempDir/expectedResults/samples/sample4/var.flt.vcf"

    assertIdenticalFiles "$tempDir/samples/sample1/consensus.fasta" "$tempDir/expectedResults/samples/sample1/consensus.fasta"
    assertIdenticalFiles "$tempDir/samples/sample2/consensus.fasta" "$tempDir/expectedResults/samples/sample2/consensus.fasta"
    assertIdenticalFiles "$tempDir/samples/sample3/consensus.fasta" "$tempDir/expectedResults/samples/sample3/consensus.fasta"
    assertIdenticalFiles "$tempDir/samples/sample4/consensus.fasta" "$tempDir/expectedResults/samples/sample4/consensus.fasta"

    assertIdenticalFiles "$tempDir/samples/sample1/reads.all.pileup" "$tempDir/expectedResults/samples/sample1/reads.all.pileup"
    assertIdenticalFiles "$tempDir/samples/sample2/reads.all.pileup" "$tempDir/expectedResults/samples/sample2/reads.all.pileup"
    assertIdenticalFiles "$tempDir/samples/sample3/reads.all.pileup" "$tempDir/expectedResults/samples/sample3/reads.all.pileup"
    assertIdenticalFiles "$tempDir/samples/sample4/reads.all.pileup" "$tempDir/expectedResults/samples/sample4/reads.all.pileup"

    assertIdenticalFiles "$tempDir/samples/sample1/consensus.vcf" "$tempDir/expectedResults/samples/sample1/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample2/consensus.vcf" "$tempDir/expectedResults/samples/sample2/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample3/consensus.vcf" "$tempDir/expectedResults/samples/sample3/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample4/consensus.vcf" "$tempDir/expectedResults/samples/sample4/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source

    assertIdenticalFiles "$tempDir/samples/sample1/metrics" "$tempDir/expectedResults/samples/sample1/metrics"
    assertIdenticalFiles "$tempDir/samples/sample2/metrics" "$tempDir/expectedResults/samples/sample2/metrics"
    assertIdenticalFiles "$tempDir/samples/sample3/metrics" "$tempDir/expectedResults/samples/sample3/metrics"
    assertIdenticalFiles "$tempDir/samples/sample4/metrics" "$tempDir/expectedResults/samples/sample4/metrics"

    assertIdenticalFiles "$tempDir/snplist_preserved.txt"                   "$tempDir/expectedResults/snplist_preserved.txt"
    assertIdenticalFiles "$tempDir/snpma_preserved.fasta"                   "$tempDir/expectedResults/snpma_preserved.fasta"
    assertIdenticalFiles "$tempDir/snpma_preserved.vcf"                     "$tempDir/expectedResults/snpma_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source --ignore-matching-lines=##bcftools
    assertIdenticalFiles "$tempDir/referenceSNP_preserved.fasta"            "$tempDir/expectedResults/referenceSNP_preserved.fasta"
    assertIdenticalFiles "$tempDir/snp_distance_pairwise_preserved.tsv"     "$tempDir/expectedResults/snp_distance_pairwise_preserved.tsv"
    assertIdenticalFiles "$tempDir/snp_distance_matrix_preserved.tsv"       "$tempDir/expectedResults/snp_distance_matrix_preserved.tsv"
    assertIdenticalFiles "$tempDir/samples/sample1/consensus_preserved.vcf" "$tempDir/expectedResults/samples/sample1/consensus_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample2/consensus_preserved.vcf" "$tempDir/expectedResults/samples/sample2/consensus_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample3/consensus_preserved.vcf" "$tempDir/expectedResults/samples/sample3/consensus_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample4/consensus_preserved.vcf" "$tempDir/expectedResults/samples/sample4/consensus_preserved.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    verifyNonEmptyReadableFile "$logDir/indexRef.log"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-1"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-2"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-4"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-3"
    verifyNonEmptyReadableFile "$logDir/callSites.log-1"
    verifyNonEmptyReadableFile "$logDir/callSites.log-2"
    verifyNonEmptyReadableFile "$logDir/callSites.log-4"
    verifyNonEmptyReadableFile "$logDir/callSites.log-3"
    verifyNonEmptyReadableFile "$logDir/mergeSites.log"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-1"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-2"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-4"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-3"
    verifyNonEmptyReadableFile "$logDir/snpMatrix.log"
    verifyNonEmptyReadableFile "$logDir/snpReference.log"
    verifyNonEmptyReadableFile "$logDir/mergeVcfs.log"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-1"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-2"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-4"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-3"
    verifyNonEmptyReadableFile "$logDir/combineMetrics.log"
    verifyNonEmptyReadableFile "$logDir/distance.log"
}


# Verify run_snp_pipeline generates correct results for the lambda data set
testRunSnpPipelineLambdaUnpaired()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Delete the 2nd half of all the paired fastq files
    rm "$tempDir"/originalInputs/samples/sample*/sample*_2.fastq

    # Simulate left-over errors from a prior run
    echo "Old errors from prior run" > "$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify correct mirroring
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample2/sample2_1.fastq" "$tempDir/samples/sample2/sample2_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample3/sample3_1.fastq" "$tempDir/samples/sample3/sample3_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample4/sample4_1.fastq" "$tempDir/samples/sample4/sample4_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/reference/lambda_virus.fasta"    "$tempDir/reference/lambda_virus.fasta"

    # Verify output results exist
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/mapReads.log-2" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/mapReads.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/mapReads.log-4" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/callSites.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/callSites.log-2" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/callSites.log-3" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/callSites.log-4" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-3" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-4" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectMetrics.log-2" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectMetrics.log-4" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectMetrics.log-3" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/mergeSites_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-3" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-4" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/distance_preserved.log" "cfsan_snp_pipeline distance finished"

}


# Verify run_snp_pipeline runs to completion with a single sample
testRunSnpPipelineLambdaSingleSample()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Delete all but one sample
    rm -rf "$tempDir"/originalInputs/samples/sample2
    rm -rf "$tempDir"/originalInputs/samples/sample3
    rm -rf "$tempDir"/originalInputs/samples/sample4

    # Simulate left-over errors from a prior run
    echo "Old errors from prior run" > "$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify correct mirroring
    assertIdenticalFiles "$tempDir/originalInputs/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_1.fastq"
    assertIdenticalFiles "$tempDir/originalInputs/reference/lambda_virus.fasta"    "$tempDir/reference/lambda_virus.fasta"

    # Verify output results exist
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    assertFileContains "$tempDir/snp_distance_pairwise.tsv" "^sample1.*sample1.*0$"
    assertFileContains "$tempDir/snp_distance_matrix.tsv" "^sample1.*0$"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"
    assertFileContains "$tempDir/snp_distance_pairwise_preserved.tsv" "^sample1.*sample1.*0$"
    assertFileContains "$tempDir/snp_distance_matrix_preserved.tsv" "^sample1.*0$"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/callSites.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/mergeSites_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/distance_preserved.log" "cfsan_snp_pipeline distance finished"

    # Verify correct results
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snpma.vcf"                     "$tempDir/samples/sample1/consensus.vcf"  # Just copy the sample VCF to the snpma.vcf
    assertIdenticalFiles "$tempDir/snpma_preserved.vcf"           "$tempDir/samples/sample1/consensus_preserved.vcf"  # Just copy the sample VCF to the snpma.vcf
}


# Verify run_snp_pipeline runs to completion with no snps and works properly when re-run
testRunSnpPipelineZeroSnps()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Clear old logs
    rm -rf $tempDir/logs*

    # Create an empty snp list and re-run the pipeline
    touch -d '-11 day' $tempDir/reference/*.fasta
    touch -d '-10 day' $tempDir/reference/*.bt2
    touch -d  '-9 day' $tempDir/samples/*/*.fastq
    touch -d  '-8 day' $tempDir/samples/*/reads.sam
    touch -d  '-7 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-6 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-4 day' $tempDir/samples/*/reads.sorted.deduped.bai
    touch -d  '-3 day' $tempDir/samples/*/realign.target.intervals
    touch -d  '-2 day' $tempDir/samples/*/reads.sorted.deduped.indelrealigned.bam
    touch -d  '-1 day' $tempDir/samples/*/reads.all.pileup
    sed -i '/PASS/d' $tempDir/samples/*/var.flt.vcf
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify output results
    verifyEmptyFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify output results exist, and no snps were found
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    lineCount=$(grep -v ">" "$tempDir/snpma.fasta" | wc -l)
    assertEquals "$tempDir/snpma.fasta should not contain any strings of bases." 0 $lineCount
    lineCount=$(grep -v ">" "$tempDir/referenceSNP.fasta" | wc -l)
    assertEquals "$tempDir/referenceSNP.fasta should not contain any strings of bases." 0 $lineCount

    verifyEmptyFile "$tempDir/snplist_preserved.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"
    lineCount=$(grep -v ">" "$tempDir/snpma_preserved.fasta" | wc -l)
    assertEquals "$tempDir/snpma_preserved.fasta should not contain any strings of bases." 0 $lineCount
    lineCount=$(grep -v ">" "$tempDir/referenceSNP_preserved.fasta" | wc -l)
    assertEquals "$tempDir/referenceSNP_preserved.fasta should not contain any strings of bases." 0 $lineCount

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/callSites.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/mergeSites_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/distance_preserved.log" "cfsan_snp_pipeline distance finished"
}


# Verify run_snp_pipeline rebuilds the snplist when at least one var.flt.vcf is missing and
# at least one var.flt.vcf is newer
testRunSnpPipelineRerunMissingVCF()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Clear old logs
    rm -rf $tempDir/logs*

    # Delete a VCF and re-run the pipeline
    touch -d '-10 day' $tempDir/reference/*.fasta
    touch -d  '-9 day' $tempDir/reference/*.bt2
    touch -d  '-8 day' $tempDir/samples/*/*.fastq
    touch -d  '-7 day' $tempDir/samples/*/reads.sam
    touch -d  '-6 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-4 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-3 day' $tempDir/samples/*/reads.sorted.deduped.bai
    touch -d  '-2 day' $tempDir/samples/*/reads.sorted.deduped.indelrealigned.bam
    touch -d  '-1 day' $tempDir/samples/*/reads.all.pileup
    touch $tempDir/samples/*/var.flt.vcf
    rm -rf "$tempDir/samples/sample1"
    sleep 1

    cfsan_snp_pipeline data configurationFile $tempDir
    echo "StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -S "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify output results
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    assertNewerFile "$tempDir/snplist.txt" "$tempDir/samples/sample2/var.flt.vcf"
    assertFileNotContains "$tempDir/snplist.txt" "sample1"

    verifyNonEmptyReadableFile "$tempDir/snplist_preserved.txt"
    assertNewerFile "$tempDir/snplist_preserved.txt" "$tempDir/samples/sample2/var.flt_preserved.vcf"
    assertFileNotContains "$tempDir/snplist_preserved.txt" "sample1"

    # Verify output results exist
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"

    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma_preserved.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP_preserved.fasta"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/callSites.log-2" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectMetrics.log-2" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/mergeSites_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/distance_preserved.log" "cfsan_snp_pipeline distance finished"
}


# Verify processing steps are skipped when output files are already fresh.
testAlreadyFreshOutputs()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> /dev/null

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Force timestamps to change so outputs are newer than inputs.
    # The files are small, quickly processed, and timestamps might not differ when we expect they will differ.
    touch -d '-15 day' $tempDir/reference/*.fasta
    touch -d '-14 day' $tempDir/reference/*.bt2
    touch -d '-13 day' $tempDir/samples/*/*.fastq
    touch -d '-12 day' $tempDir/samples/*/reads.sam
    touch -d '-11 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d '-10 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-9 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-8 day' $tempDir/samples/*/reads.sorted.deduped.bai
    touch -d  '-7 day' $tempDir/samples/*/realign.target.intervals
    touch -d  '-6 day' $tempDir/samples/*/reads.sorted.deduped.indelrealigned.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-4 day' $tempDir/samples/*/var.flt.vcf
    touch -d  '-3 day' $tempDir/samples/*/var.flt_preserved.vcf
    touch -d  '-3 day' $tempDir/samples/*/var.flt_removed.vcf
    touch -d  '-2 day' $tempDir/snplist.txt
    touch -d  '-2 day' $tempDir/snplist_preserved.txt
    touch -d  '-1 day' $tempDir/samples/*/consensus*.vcf

    # Test special cfsan_snp_pipeline collect_metrics result persistence
    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=20000"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=110"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=94.55"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=286.84"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=23.21"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=46"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=32"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=46"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=32"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=0"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=0"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample="
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSamplePreserved="

    echo numberReads=AA > "$tempDir/samples/sample1/metrics"
    echo numberDupReads=BB >> "$tempDir/samples/sample1/metrics"
    echo percentReadsMapped=CC >> "$tempDir/samples/sample1/metrics"
    echo aveInsertSize=DD >> "$tempDir/samples/sample1/metrics"
    echo avePileupDepth=EE >> "$tempDir/samples/sample1/metrics"
    echo phase1Snps=FF >> "$tempDir/samples/sample1/metrics"
    echo phase1SnpsPreserved=GG >> "$tempDir/samples/sample1/metrics"
    echo snps=HH >> "$tempDir/samples/sample1/metrics"
    echo snpsPreserved=II >> "$tempDir/samples/sample1/metrics"
    echo missingPos=JJ >> "$tempDir/samples/sample1/metrics"
    echo missingPosPreserved=KK >> "$tempDir/samples/sample1/metrics"

    # Remove unwanted log files
    rm -rf $tempDir/logs*

    # Re-run the pipeline
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> /dev/null

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/indexRef.log"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-1"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-2"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-4"
    verifyNonEmptyReadableFile "$logDir/mapReads.log-3"
    verifyNonEmptyReadableFile "$logDir/callSites.log-1"
    verifyNonEmptyReadableFile "$logDir/callSites.log-2"
    verifyNonEmptyReadableFile "$logDir/callSites.log-4"
    verifyNonEmptyReadableFile "$logDir/callSites.log-3"
    verifyNonEmptyReadableFile "$logDir/mergeSites.log"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-1"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-2"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-4"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-3"
    verifyNonEmptyReadableFile "$logDir/snpMatrix.log"
    verifyNonEmptyReadableFile "$logDir/snpReference.log"
    verifyNonEmptyReadableFile "$logDir/mergeVcfs.log"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-1"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-2"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-4"
    verifyNonEmptyReadableFile "$logDir/collectMetrics.log-3"

    assertFileContains "$logDir/indexRef.log" "lambda_virus.rev.1.bt2 is already freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/indexRef.log" "lambda_virus.fasta.fai is already freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/indexRef.log" "lambda_virus.dict is already freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mapReads.log-1" "sample1 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "sample2 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "sample4 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "sample3 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Unsorted bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Unsorted bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Unsorted bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Unsorted bam file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Sorted bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Sorted bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Sorted bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Sorted bam file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Deduped bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Deduped bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Deduped bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Deduped bam file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Bam file index is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Bam file index is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Bam file index is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Bam file index is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Realign targets file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Realign targets file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Realign targets file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Realign targets file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-1" "Indelrealigned bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-2" "Indelrealigned bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-3" "Indelrealigned bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/mapReads.log-4" "Indelrealigned bam file is already freshly created for sample3.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/callSites.log-1" "Pileup file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-2" "Pileup file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-3" "Pileup file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-4" "Pileup file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-1" "VCF file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-2" "VCF file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-3" "VCF file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callSites.log-4" "VCF file is already freshly created for sample3.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeSites.log" "snplist.txt has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/callConsensus.log-1" "sample1/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-2" "sample2/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-3" "sample4/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-4" "sample3/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpMatrix.log" "/snpma.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpReference.log" "referenceSNP.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeVcfs.log" "Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."

    assertFileNotContains "$logDir/collectMetrics.log-1" "already freshly created"
    assertFileNotContains "$logDir/collectMetrics.log-2" "already freshly created"
    assertFileNotContains "$logDir/collectMetrics.log-3" "already freshly created"
    assertFileNotContains "$logDir/collectMetrics.log-4" "already freshly created"

    assertFileContains "$logDir/distance.log" "have already been freshly built.  Use the -f option to force a rebuild"

    # =======
    assertFileContains "$logDir/filterRegions.log" "All preserved and removed vcf files are already freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeSites_preserved.log" "snplist_preserved.txt has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/callConsensus_preserved.log-1" "sample1/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-2" "sample2/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-3" "sample4/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-4" "sample3/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpMatrix_preserved.log" "/snpma_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpReference_preserved.log" "referenceSNP_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeVcfs_preserved.log" "Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/distance_preserved.log" "have already been freshly built.  Use the -f option to force a rebuild"

    # Special cfsan_snp_pipeline collect_metrics re-use last metrics
    assertFileNotContains "$tempDir/samples/sample1/metrics" "numberReads=20000"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "numberDupReads=110"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=94.55"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "aveInsertSize=286.84"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "avePileupDepth=23.22"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "phase1Snps=46"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=32"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "snps=46"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "snpsPreserved=32"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "missingPos=0"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=0"

    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=AA"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=BB"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=CC"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=DD"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=EE"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=FF"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=GG"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=HH"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=II"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=JJ"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=KK"

    assertFileContains "$tempDir/metrics.tsv" "sample1.*AA.*BB.*CC.*DD.*EE.*FF.*GG.*HH.*II.*JJ.*KK"
}


# Verify underscores in metrics.tsv column headers can be controlled with configuration file
testRunSnpPipelineMetricsColumnHeadingsUnderscores()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify output metrics
    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/metrics.tsv" "$tempDir/expectedResults/metrics.tsv"
    head -n 1 "$tempDir/metrics.tsv" | grep "_" > /dev/null
    assertTrue "No underscores were found in the metrics column headings"  $?

    # Delete the metrics file and re-run with the option to use spaces
    rm "$tempDir/metrics.tsv"
    cfsan_snp_pipeline data configurationFile $tempDir
    echo 'CombineMetrics_ExtraParams="-s"' >> "$tempDir/snppipeline.conf"
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify output metrics have no underscores
    head -n 1 "$tempDir/metrics.tsv" | grep "_" > /dev/null
    assertFalse "Underscores should not be found in the metrics column headings when using -s combine_metrics option"  $?
}


# Verify samples with excessive snps are excluded from the snplist, snp matrix, and snpma.vcf.
testRunSnpPipelineExcessiveSnps()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Create a config file with a low enough maxsnps setting to block a sample
    cfsan_snp_pipeline data configurationFile $tempDir
    sed -i s:MaxSnps=-1:MaxSnps=40: "$tempDir/snppipeline.conf"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

	# Verify each pipeline stage runs to completion
    logDir=$(echo $(ls -d $tempDir/logs*))
    assertFileContains "$logDir/indexRef.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/mapReads.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/callSites.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/mergeSites.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/distance.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterRegions.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/mergeSites_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcfs_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/distance_preserved.log" "cfsan_snp_pipeline distance finished"

    # Verify output
    # After removing the abnormal high-density snps, sample1 has more than 40 snps, so it is included in the analysis - non-preserved only
    # After removing the abnormal high-density snps, sample2 has more than 40 snps, so it is included in the analysis - both non-preserved and preserved
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample=Excluded$"
    assertFileContains "$tempDir/samples/sample2/metrics" "excludedSample=Excluded$"
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSamplePreserved=$"
    assertFileContains "$tempDir/samples/sample2/metrics" "excludedSamplePreserved=Excluded$"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1Snps=46"
    assertFileContains "$tempDir/samples/sample2/metrics" "phase1Snps=44"
    assertFileContains "$tempDir/samples/sample1/metrics" "phase1SnpsPreserved=32"
    assertFileContains "$tempDir/samples/sample2/metrics" "phase1SnpsPreserved=41"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=$"
    assertFileContains "$tempDir/samples/sample2/metrics" "snps=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snpsPreserved=32$"
    assertFileContains "$tempDir/samples/sample2/metrics" "snpsPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample2/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPosPreserved=0$"
    assertFileContains "$tempDir/samples/sample2/metrics" "missingPosPreserved=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "errorList=.*Excluded: exceeded 40 maxsnps"
    assertFileContains "$tempDir/samples/sample2/metrics" "errorList=.*Excluded: exceeded 40 maxsnps.*Excluded: preserved exceeded 40 maxsnps"
    assertFileContains "$tempDir/metrics.tsv"             "sample1.*Excluded.*Excluded: exceeded 40 maxsnps"
    assertFileContains "$tempDir/metrics.tsv"             "sample2.*Excluded.*Excluded: exceeded 40 maxsnps.*Excluded: preserved exceeded 40 maxsnps"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "Consensus.*not found"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "Consensus.*not found"
    assertFileNotContains "$tempDir/metrics.tsv"             "sample1.*Consensus.*not found"

    assertFileNotContains "$tempDir/snplist.txt" "sample1"
    assertFileNotContains "$tempDir/snpma.fasta" "sample1"
    assertFileNotContains "$tempDir/snpma.vcf"   "sample1"
    assertFileNotContains "$tempDir/snp_distance_pairwise.tsv" "sample1"
    assertFileNotContains "$tempDir/snp_distance_matrix.tsv" "sample1"

    assertFileNotContains "$tempDir/snplist.txt" "sample2"
    assertFileNotContains "$tempDir/snpma.fasta" "sample2"
    assertFileNotContains "$tempDir/snpma.vcf"   "sample2"
    assertFileNotContains "$tempDir/snp_distance_pairwise.tsv" "sample2"
    assertFileNotContains "$tempDir/snp_distance_matrix.tsv" "sample2"

    assertFileContains "$tempDir/snplist_preserved.txt" "sample1"
    assertFileContains "$tempDir/snpma_preserved.fasta" "sample1"
    assertFileContains "$tempDir/snpma_preserved.vcf"   "sample1"
    assertFileContains "$tempDir/snp_distance_pairwise_preserved.tsv" "sample1"
    assertFileContains "$tempDir/snp_distance_matrix_preserved.tsv" "sample1"

    assertFileNotContains "$tempDir/snplist_preserved.txt" "sample2"
    assertFileNotContains "$tempDir/snpma_preserved.fasta" "sample2"
    assertFileNotContains "$tempDir/snpma_preserved.vcf"   "sample2"
    assertFileNotContains "$tempDir/snp_distance_pairwise_preserved.tsv" "sample2"
    assertFileNotContains "$tempDir/snp_distance_matrix_preserved.tsv" "sample2"

    cfsan_snp_pipeline data lambdaVirusExpectedResults $tempDir/expectedResults
    grep -v sample[12] "$tempDir/expectedResults/metrics.tsv" > "$tempDir/expectedResults/metrics.withoutSample1and2.tsv"
    grep -v sample[12] "$tempDir/metrics.tsv" > "$tempDir/metrics.withoutSample1and2.tsv"
    assertIdenticalFiles "$tempDir/metrics.withoutSample1and2.tsv" "$tempDir/expectedResults/metrics.withoutSample1and2.tsv"
}


# Verify run_snp_pipeline generates correct results for the Salmonella Agona data set
testRunSnpPipelineAgona()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    cfsan_snp_pipeline data agonaInputs "$tempDir"

    # Download sample data from SRA at NCBI.
     mkdir "$tempDir/samples"
    < "$tempDir/sampleList" xargs -P 10 -I % sh -c "mkdir $tempDir/samples/%; fastq-dump --gzip --origfmt --split-files --outdir $tempDir/samples/% % > /dev/null;"

    # Verify the agona data was downloaded properly
    cd "$tempDir"
    sha256sum -c sha256sumCheck &> /dev/null
    rc=$?
    assertTrue "The agona data set was not downloaded properly" "$rc"
    cd - &> /dev/null
    if [[ $rc != 0 ]]; then return $rc; fi  # skip the long-running agona pipeline run

    # Simulate left-over errors from a prior run
    echo "Old errors from prior run" > "$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference//NC_011149.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "There were errors processing some samples"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify correct results
    cfsan_snp_pipeline data agonaExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snplist.txt"        "$tempDir/expectedResults/snplist.txt"
    assertIdenticalFiles "$tempDir/snpma.fasta"        "$tempDir/expectedResults/snpma.fasta"
    assertIdenticalFiles "$tempDir/snpma.vcf"          "$tempDir/expectedResults/snpma.vcf" "--ignore-matching-lines=##fileDate" "--ignore-matching-lines=##source" "--ignore-matching-lines=##bcftools"
    assertIdenticalFiles "$tempDir/referenceSNP.fasta" "$tempDir/expectedResults/referenceSNP.fasta"

    assertIdenticalFiles "$tempDir/metrics.tsv"               "$tempDir/expectedResults/metrics.tsv"
    assertIdenticalFiles "$tempDir/snp_distance_pairwise.tsv" "$tempDir/expectedResults/snp_distance_pairwise.tsv"
    assertIdenticalFiles "$tempDir/snp_distance_matrix.tsv"   "$tempDir/expectedResults/snp_distance_matrix.tsv"

    assertIdenticalFiles "$tempDir/snplist_preserved.txt"             "$tempDir/expectedResults/snplist_preserved.txt"
    assertIdenticalFiles "$tempDir/snpma_preserved.fasta"             "$tempDir/expectedResults/snpma_preserved.fasta"
    assertIdenticalFiles "$tempDir/snpma_preserved.vcf"               "$tempDir/expectedResults/snpma_preserved.vcf" "--ignore-matching-lines=##fileDate" "--ignore-matching-lines=##source" "--ignore-matching-lines=##bcftools"
    assertIdenticalFiles "$tempDir/referenceSNP_preserved.fasta"      "$tempDir/expectedResults/referenceSNP_preserved.fasta"
    assertIdenticalFiles "$tempDir/snp_distance_matrix_preserved.tsv" "$tempDir/expectedResults/snp_distance_matrix_preserved.tsv"

    # Verify proper read group tags
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "@RG"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "ID:H92ARADXX.2"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "SM:SRR1566386"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "LB:1"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "PL:illumina"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "PU:H92ARADXX.2.SRR1566386"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "RG:Z:H92ARADXX.2"
    assertFileContains "$tempDir/samples/SRR1566386/reads.sam" "HWI-ST1029:112:H92ARADXX:2:1101:1623:2163"

    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "@RG"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "ID:A87H3.1"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "SM:SRR2566901"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "LB:1"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "PL:illumina"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "PU:A87H3.1.SRR2566901"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "RG:Z:A87H3.1"
    assertFileContains "$tempDir/samples/SRR2566901/reads.sam" "M02070:11:000000000-A87H3:1:1101:16393:1371"
}


# load shunit2 and execute all the tests in this script
. shunit2
