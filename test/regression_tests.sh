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
    [ "$1" -nt "$2" ]
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
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "prepReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "alignSampleToReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "prepSamples.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "snp_filter.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_list.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "call_consensus.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_matrix.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_reference_seq.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "collectSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "combineSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "calculate_snp_distances.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "$aligner"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "samtools"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "java"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "tabix"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bgzip"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bcftools"
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to VarScan"
    assertFileContains "$tempDir/error.log" "Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html"

    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "cfsan_snp_pipeline"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "prepReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "alignSampleToReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "prepSamples.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "snp_filter.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_list.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "call_consensus.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_matrix.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_reference_seq.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "collectSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "combineSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "calculate_snp_distances.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "$aligner"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "samtools"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "java"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "tabix"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "bgzip"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "bcftools"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to VarScan"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html"
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyPicardRequiredRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_RemoveDuplicateReads=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to Picard"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to Picard"
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyPicardNotRequiredRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_RemoveDuplicateReads=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
    assertFileNotContains "$tempDir/error.log" "CLASSPATH is not configured with the path to Picard"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "CLASSPATH is not configured with the path to Picard"
}

# Verify the index_ref command detects a misconfigured environment variable
tryPrepReferenceEnvironmentRaiseGlobalError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref failed"
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/prepReference.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
tryPrepReferenceEmptyFastaFileRaiseGlobalError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the reference file was empty." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref failed"
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileContains "$logDir/prepReference.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}

# Verify the index_ref command detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}


# Verify the index_ref command detects bowtie error and emits the global error marker file.
tryPrepReferenceBowtieIndexTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref / bowtie returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceBowtieIndexTrap 100
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceBowtieIndexTrap 100
}

# Verify the index_ref command detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceBowtieIndexTrap 100
}


# Verify the index_ref command detects smalt error and emits the global error marker file.
tryPrepReferenceSmaltIndexTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline index_ref / smalt returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/error.log" "smalt index"
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the index_ref command detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
tryPrepReferenceSamtoolsFaidxTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
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

    # SAMtools 0.1.19 prepReference.log
    # ---------------------------------
    # [fai_build_core] different line length in sequence 'gi|9626243|ref|NC_001416.1|'.
    # /home/steven.davis/.virtualenvs/snp-pipeline-2.7/bin/prepReference.sh: line 176: 30719 Segmentation fault      (core dumped) samtools faidx $SamtoolsFaidx_ExtraParams "$referenceFilePath"

    # SAMtools 1.3 prepReference.log
    # --------------------------
    # [fai_build_core] different line length in sequence 'gi|9626243|ref|NC_001416.1|'.
    # Error: /tmp/shunit.VB9zhq/tmp/tmp.qQsT4ORbSy/reference/lambda_virus.fasta.fai does not exist after running samtools faidx.

    assertEquals "cfsan_snp_pipeline index_ref returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline index_ref"
    assertFileContains "$tempDir/error.log" "samtools faidx"
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileNotContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceSamtoolsFaidxTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceSamtoolsFaidxTrap 100
}

# Verify the index_ref command detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceSamtoolsFaidxTrap 100
}



# Verify the map_reads command detects a misconfigured environment variable.
tryAlignSampleToReferenceEnvironmentRaiseGlobalError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately misconfigure the environment
    export SnpPipeline_Aligner=garbage

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/alignSampleToReference.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the map_reads command detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}

# Verify the map_reads command detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}

# Verify the map_reads command detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}


# Verify the map_reads command detects an Missing reference file
tryAlignSampleToReferenceMissingReferenceRaiseGlobalError()
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
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"
}

# Verify the map_reads command detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the map_reads command detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the map_reads command detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}


# Verify the map_reads command detects a missing sample file
tryAlignSampleToReferenceMissingSample1RaiseSampleError()
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
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 100
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 98
}

# Verify the map_reads command detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 100
}


# Verify the map_reads command detects a missing sample file
tryAlignSampleToReferenceMissingSample2RaiseSampleError()
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
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline map_reads failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 100
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 98
}

# Verify the map_reads command detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 100
}


# Verify the map_reads command detects bowtie alignment error.
tryAlignSampleToReferenceBowtieAlignTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log-1"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/alignSampleToReference.log-1" "bowtie2"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/alignSampleToReference.log-3"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/alignSampleToReference.log-3" "bowtie2"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Use the -f option to force a rebuild"
}


# Verify the map_reads command detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceBowtieAlignTrap 100
}

# Verify the map_reads command detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceBowtieAlignTrap 98
}

# Verify the map_reads command detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceBowtieAlignTrap 100
}


# Verify the map_reads command detects smalt alignment error.
tryAlignSampleToReferenceSmaltAlignTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log-1"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/alignSampleToReference.log-1" "smalt map"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/alignSampleToReference.log-3"
    errorCode=$?

    # Verify the map_reads command error handling behavior
    assertEquals "cfsan_snp_pipeline map_reads with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Error detected while running cfsan_snp_pipeline map_reads."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/alignSampleToReference.log-3" "smalt map"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}


# Verify the map_reads command detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceSmaltAlignTrap 100
}

# Verify the map_reads command detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceSmaltAlignTrap 98
}

# Verify the map_reads command detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceSmaltAlignTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects a Missing reference file
tryPrepSamplesMissingReferenceRaiseGlobalError()
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
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta" "xxxx" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/prepSamples.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing Reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a missing reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline call_sites script detects a missing sample sam file
tryPrepSamplesMissingSamFileRaiseSampleError()
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
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when the reads.sam file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam does not exist"
    assertFileContains "$logDir/prepSamples.log" "Sample SAM file $tempDir/samples/sample1/reads.sam does not exist"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the cfsan_snp_pipeline call_sites script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesMissingSamFileRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesMissingSamFileRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline call_sites script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesMissingSamFileRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline call_sites script detects Samtools view failure.
tryPrepSamplesSamtoolsViewTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the SAM file and re-run cfsan_snp_pipeline call_sites
    echo "Garbage" > "$tempDir/samples/sample1/reads.sam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when the input SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileContains "$tempDir/error.log" "samtools view"
    assertFileContains "$logDir/prepSamples.log" "samtools view"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsViewTrap 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsViewTrap 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsViewTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects Samtools sort failure.
tryPrepSamplesSamtoolsSortTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the reads.unsorted.bam file and re-run cfsan_snp_pipeline call_sites
    echo "Garbage" > "$tempDir/samples/sample1/reads.unsorted.bam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when reads.unsorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileContains "$tempDir/error.log" "samtools sort"
    assertFileContains "$logDir/prepSamples.log" "samtools sort"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Sorted bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsSortTrap 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsSortTrap 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsSortTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects Picard MarkDuplicates failure.
tryPrepSamplesPicardMarkDuplicatesTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the reads.sorted.bam file and re-run cfsan_snp_pipeline call_sites
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.bam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when reads.sorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileContains "$tempDir/error.log" "picard"
    assertFileContains "$tempDir/error.log" "MarkDuplicates"
    assertFileContains "$logDir/prepSamples.log" "picard"
    assertFileContains "$logDir/prepSamples.log" "MarkDuplicates"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Pileup file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Picard MarkDuplicates failure.
testPrepSamplesPicardMarkDuplicatesTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesPicardMarkDuplicatesTrap 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Picard MarkDuplicates failure.
testPrepSamplesPicardMarkDuplicatesTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesPicardMarkDuplicatesTrap 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Picard MarkDuplicates failure.
testPrepSamplesPicardMarkDuplicatesTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesPicardMarkDuplicatesTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Picard MarkDuplicates.
tryPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Clear CLASSPATH and re-run cfsan_snp_pipeline call_sites
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when CLASSPATH unset failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute Picard. Define the path to Picard in the CLASSPATH environment variable."
    assertFileContains "$logDir/prepSamples.log" "Error: cannot execute Picard. Define the path to Picard in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Deduped bam file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run PicardMarkDuplicates.
testPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run PicardMarkDuplicates.
testPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run PicardMarkDuplicates.
testPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesPicardMarkDuplicatesClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
tryPrepSamplesSamtoolsMpileupTrap()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the reads.sorted.deduped.bam file and re-run cfsan_snp_pipeline call_sites
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when reads.sorted.deduped.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running cfsan_snp_pipeline call_sites."
    assertFileContains "$tempDir/error.log" "samtools mpileup"
    assertFileContains "$logDir/prepSamples.log" "samtools mpileup"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Pileup file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsMpileupTrap 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsMpileupTrap 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsMpileupTrap 100
}


# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
tryPrepSamplesVarscanRaiseSampleError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Configure VarScan with invalid parameter settings and re-run cfsan_snp_pipeline call_sites
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.all.pileup"
    export VarscanMpileup2snp_ExtraParams="--min-coverage -1 --min-reads 99999999 --min-avg_qual -100 --min-var-freq 2 --output-vcf 2 --invalid-parameter"
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?
    export VarscanMpileup2snp_ExtraParams=""

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when varscan failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/prepSamples.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/prepSamples.log" "VarScan mpileup2snp"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Vcf file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesVarscanRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesVarscanRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline call_sites script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesVarscanRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
tryPrepSamplesVarscanClasspathRaiseGlobalError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run cfsan_snp_pipeline call_sites normally
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Configure VarScan with invalid parameter settings and re-run cfsan_snp_pipeline call_sites
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.deduped.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.all.pileup"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    rm "$tempDir/error.log" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify cfsan_snp_pipeline call_sites error handling behavior
    assertEquals "cfsan_snp_pipeline call_sites returned incorrect error code when CLASSPATH unset failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_sites failed"
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileContains "$logDir/prepSamples.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/prepSamples.log" "cfsan_snp_pipeline call_sites finished"
    assertFileNotContains "$logDir/prepSamples.log" "Vcf file is already freshly created"
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the cfsan_snp_pipeline call_sites script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
trySnpFilterPermissionTrap()
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
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_preserved.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for preserved SNPs"
    assertFileNotContains "$logDir/snp_filter.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_preserved.vcf"
    rm "$tempDir/error.log"

    # 2) var.flt_removed.vcf =========
    # Make the output file var.flt_removed.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_removed.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_removed.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_removed.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_removed.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for removed SNPs"
    assertFileNotContains "$logDir/snp_filter.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_removed.vcf"
    rm "$tempDir/error.log"
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testSnpFilterPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterPermissionTrap 100
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testSnpFilterPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false

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
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_preserved.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for preserved SNPs"
    assertFileNotContains "$logDir/snp_filter.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_preserved.vcf"
    rm "$tempDir/error.log"

    # 2) var.flt_removed.vcf =========
    # Make the output file var.flt_removed.vcf unwritable
    mkdir -p "$tempDir/samples/sample1"
    touch "$tempDir/samples/sample1/var.flt_removed.vcf"
    chmod -w "$tempDir/samples/sample1/var.flt_removed.vcf"

    # Try to run cfsan_snp_pipeline filter_regions -- it should have problems writing to sample1/var.flt_removed.vcf
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline filter_regions error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when var.flt_removed.vcf was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Cannot create the file for removed SNPs"
    assertFileNotContains "$logDir/snp_filter.log" "Error detected while running cfsan_snp_pipeline filter_regions"
    assertFileContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/samples/sample1/var.flt_removed.vcf"
    rm "$tempDir/error.log"
}

# Verify the cfsan_snp_pipeline filter_regions script traps attempts to write to unwritable file
testSnpFilterPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterPermissionTrap 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
trySnpFilterMissingSampleDirRaiseGlobalError()
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
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/snp_filter.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testSnpFilterMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testSnpFilterMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    trySnpFilterMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing sample directories file
testSnpFilterMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
trySnpFilterMissingReferenceRaiseGlobalError()
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
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/non-exist-reference" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/non-exist-reference does not exist"
    assertFileContains "$logDir/snp_filter.log" "Reference file $tempDir/non-exist-reference does not exist"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testSnpFilterMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testSnpFilterMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    trySnpFilterMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing reference file
testSnpFilterMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
trySnpFilterMissingOutgroupRaiseGlobalError()
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
    cfsan_snp_pipeline filter_regions -g "$tempDir/outgroup" "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when file of outgroup samples file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "File of outgroup samples $tempDir/outgroup does not exist"
    assertFileContains "$logDir/snp_filter.log" "File of outgroup samples $tempDir/outgroup does not exist"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testSnpFilterMissingOutgroupRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterMissingOutgroupRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testSnpFilterMissingOutgroupRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    trySnpFilterMissingOutgroupRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects missing outgroup samples file
testSnpFilterMissingOutgroupRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterMissingOutgroupRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
trySnpFilterMissingVcfRaiseGlobalError()
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
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when all var.flt.vcf were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$logDir/snp_filter.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testSnpFilterMissingVcfRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterMissingVcfRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testSnpFilterMissingVcfRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    trySnpFilterMissingVcfRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline filter_regions script detects all vcf files missing
testSnpFilterMissingVcfRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterMissingVcfRaiseGlobalError 100
}


# Verify the filter_regions command detects missing some VCF files, but not all
trySnpFilterMissingVcfRaiseSampleError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline filter_regions -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when some var.flt.vcf were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/snp_filter.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the filter_regions command detects missing some VCF files, but not all
testSnpFilterMissingVcfRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    trySnpFilterMissingVcfRaiseSampleError 100
}

# Verify the filter_regions command detects missing some VCF files, but not all
testSnpFilterMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline filter_regions -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirList.txt" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/snp_filter.log"
    errorCode=$?

    # Verify the filter_regions command error handling behavior
    assertEquals "cfsan_snp_pipeline filter_regions returned incorrect error code when some var.flt.vcf were missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileNotContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions failed."
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/snp_filter.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/snp_filter.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/snp_filter.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileNotContains "$logDir/snp_filter.log" "Use the -f option to force a rebuild"
}

# Verify the filter_regions command detects missing some VCF files, but not all
testSnpFilterMissingVcfRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    trySnpFilterMissingVcfRaiseSampleError 100
}


# Verify the filter_regions command uses all the input vcf files to produce the outputs
# even when some of the samples are already fresh.
testSnpFilterPartialRebuild()
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
    touch -d  '-5 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-4 day' $tempDir/samples/*/var.flt.vcf
    touch -d  '-3 day' $tempDir/samples/*/var.flt_preserved.vcf
    touch -d  '-3 day' $tempDir/samples/*/var.flt_removed.vcf

    # Remove unwanted log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    rm -rf $logDir
    mkdir -p $logDir

    # Remove the results for one of the samples
    rm $tempDir/samples/sample1/var.flt_preserved.vcf

    # Re-run cfsan_snp_pipeline filter_regions -- this should only rebuild results for sample1, but it should use the var.flt.vcf input file for all samples
    cfsan_snp_pipeline filter_regions "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" > "$logDir/filterAbnormalSNP.log"

    # Verify log files
    verifyNonEmptyReadableFile "$logDir/filterAbnormalSNP.log"
    assertFileNotContains "$logDir/filterAbnormalSNP.log" "already freshly built"

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
testSnpFilterOutgroup()
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
    touch -d  '-4 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-3 day' $tempDir/samples/*/var.flt.vcf

    # Remove unwanted log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    rm -rf $logDir
    mkdir -p $logDir

    # One of the samples is an outgroup
    outgroup="sample4" # this test only works when sample 4 is the outgroup
    echo $outgroup > "$tempDir/outgroup.txt"

    # Re-run cfsan_snp_pipeline filter_regions --
    cfsan_snp_pipeline filter_regions --out_group "$tempDir/outgroup.txt" "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" > "$logDir/filterAbnormalSNP.log"

    # Verify log files
    verifyNonEmptyReadableFile "$logDir/filterAbnormalSNP.log"
    assertFileNotContains "$logDir/filterAbnormalSNP.log" "already freshly built"
    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"

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
tryCreateSnpListPermissionTrap()
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
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline merge_sites"
    assertFileNotContains "$logDir/create_snp_list.log" "Error detected while running cfsan_snp_pipeline merge_sites"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snplist.txt"
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListPermissionTrap 100
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListPermissionTrap 100
}

# Verify the cfsan_snp_pipeline merge_sites script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListPermissionTrap 100
}



# Verify the merge_sites command detects failure.
tryCreateSnpListMissingSampleDirRaiseGlobalError()
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
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/create_snp_list.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}


# Verify the merge_sites command detects failure.
tryCreateSnpListMissingVcfRaiseGlobalError()
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
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "Error: all 4 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}


# Verify the merge_sites command detects failure.
tryCreateSnpListMissingVcfRaiseSampleError()
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline merge_sites -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingVcfRaiseSampleError 100
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline merge_sites -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline merge_sites -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" "$tempDir/sampleDirectories.txt.OrigVCF.filtered" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify merge_sites error handling behavior
    assertEquals "cfsan_snp_pipeline merge_sites returned incorrect error code when var.flt.vcf was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileNotContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites failed."
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the merge_sites command detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingVcfRaiseSampleError 100
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
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline call_consensus"
    assertFileNotContains "$logDir/call_consensus.log" "Error detected while running cfsan_snp_pipeline call_consensus"
    assertFileContains "$tempDir/error.log" "function read_snp_position_list at line"
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusCorruptSnplistTrap 100
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusCorruptSnplistTrap 98
}

# Verify the call_consensus command traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
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
    cfsan_snp_pipeline call_consensus -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the snplist file"
    assertFileContains "$logDir/call_consensus.log" "cannot call consensus without the snplist file"
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/call_consensus.log" "Snplist file snplist.txt does not exist"
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus command detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
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
    cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline call_consensus with empty snplist
    touch "$tempDir/snplist.txt"
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/samples/sample1/consensus.fasta"  --vcfFileName "$tempDir/samples/sample1/consensus.vcf" "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when snplist was empty." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.vcf"
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/call_consensus.log" "cannot call consensus without the snplist file"
    assertFileNotContains "$logDir/call_consensus.log" "Snplist file $tempDir/snplist.txt is empty"
    assertFileContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus command does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStopUnset()
{
    unset SnpPipeline_StopOnSampleError
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
    cfsan_snp_pipeline call_consensus -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when pileup file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the pileup file"
    assertFileContains "$logDir/call_consensus.log" "cannot call consensus without the pileup file"
    assertFileContains "$tempDir/error.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileContains "$logDir/call_consensus.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusMissingPileupRaiseSampleError 100
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusMissingPileupRaiseSampleError 98
}

# Verify the call_consensus command detects missing pileup file.
testCallConsensusMissingPileupRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
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
    cfsan_snp_pipeline call_consensus -e "$tempDir/samples/sample1/excludeFile.vcf" -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify the call_consensus command error handling behavior
    assertEquals "cfsan_snp_pipeline call_consensus returned incorrect error code when exclude file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the file of excluded positions"
    assertFileContains "$logDir/call_consensus.log" "cannot call consensus without the file of excluded positions"
    assertFileContains "$tempDir/error.log" "Exclude file $tempDir/samples/sample1/excludeFile.vcf does not exist"
    assertFileContains "$logDir/call_consensus.log" "Exclude file $tempDir/samples/sample1/excludeFile.vcf does not exist"
    assertFileNotContains "$logDir/call_consensus.log" "cfsan_snp_pipeline call_consensus finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusMissingExcludeRaiseSampleError 100
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusMissingExcludeRaiseSampleError 98
}

# Verify the call_consensus command detects missing exclude file.
testCallConsensusMissingExcludeRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCallConsensusMissingExcludeRaiseSampleError 100
}


# Verify the merge_vcfs command detects failure.
tryMergeVcfCorruptVcfTrap()
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
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when consensus.vcf was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline merge_vcfs."
    assertFileNotContains "$logDir/mergeVcf.log" "Error detected while running cfsan_snp_pipeline merge_vcfs."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
testMergeVcfCorruptVcfTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfCorruptVcfTrap 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfCorruptVcfTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfCorruptVcfTrap 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfCorruptVcfTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfCorruptVcfTrap 100
}


# Verify the merge_vcfs command detects failure.
tryMergeVcfMissingSampleDirRaiseGlobalError()
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
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist."
    assertFileContains "$logDir/mergeVcf.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}


# Verify the merge_vcfs command detects a missing consensus VCF file.
tryMergeVcfMissingVcfRaiseSampleError()
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
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when vcf file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects a missing consensus VCF file.
testMergeVcfMissingVcfRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfMissingVcfRaiseSampleError 100
}

# Verify the merge_vcfs command detects a missing consensus VCF file - but continues running.
testMergeVcfMissingVcfRaiseSampleErrorNoStop()
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
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with missing vcf files
    rm "$tempDir/samples/sample1/consensus.vcf"
    rm "$tempDir/snpma.vcf"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command keeps running when only one vcf file is missing
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when vcf file was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
}

# Verify the merge_vcfs command detects a missing consensus VCF file.
testMergeVcfMissingVcfRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfMissingVcfRaiseSampleError 100
}


# Verify the merge_vcfs command simply copies the input consensus VCF file when there is only one sample
testMergeVcfOnlyOneSample()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=true
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline merge_vcfs with only one vcf file
    mkdir -p "$tempDir/samples/sample1"
    echo "$tempDir/samples/sample1" > "$tempDir/sampleDirectories.txt"
    echo "Dummy VCF contents" > "$tempDir/samples/sample1/consensus.vcf"
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command copies the input consensus VCF file when there is only one sample
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when there was only one vcf file." 0 $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcf.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    assertIdenticalFiles "$tempDir/samples/sample1/consensus.vcf" "$tempDir/snpma.vcf"
}


# Verify the merge_vcfs command detects all the consensus VCF files missing
tryMergeVcfZeroGoodSamplesRaiseGlobalError()
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
    cfsan_snp_pipeline merge_vcfs -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify the merge_vcfs command error handling behavior
    assertEquals "cfsan_snp_pipeline merge_vcfs returned incorrect error code when no good VCF files." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs failed."
    assertFileContains "$tempDir/error.log" "There are no vcf files to merge."
    assertFileContains "$logDir/mergeVcf.log" "There are no vcf files to merge."
    assertFileNotContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the merge_vcfs command detects failure.
#testMergeVcfZeroGoodSamplesRaiseGlobalErrorStop()
#{
#    # Nothing to test, SnpPipeline_StopOnSampleError must be false to test this code path.
#    # Otherwise, the first missing VCF file will trigger stop upon sample error -- already tested.
#}

# Verify the merge_vcfs command detects failure.
testMergeVcfZeroGoodSamplesRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfZeroGoodSamplesRaiseGlobalError 100
}

# Verify the merge_vcfs command detects failure.
#testMergeVcfZeroGoodSamplesRaiseGlobalErrorStopUnset()
#{
#    # Nothing to test, SnpPipeline_StopOnSampleError must be false to test this code path.
#    # Otherwise, the first missing VCF file will trigger stop upon sample error -- already tested.
#}



# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
tryCreateSnpMatrixPermissionTrap()
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
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline snp_matrix"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Error detected while running cfsan_snp_pipeline snp_matrix"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snpma.fasta"
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_matrix script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixPermissionTrap 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
tryCreateSnpMatrixMissingSampleDirRaiseGlobalError()
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
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
tryCreateSnpMatrixMissingConsensusRaiseGlobalError()
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
    #cfsan_snp_pipeline index_ref "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    #cfsan_snp_pipeline map_reads "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    #cfsan_snp_pipeline call_sites "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run cfsan_snp_pipeline snp_matrix -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpmap.fasta"  "$tempDir/sampleDirList.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when all consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample2/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample3/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample2/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample3/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: all 4 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "Error: all 4 consensus fasta files were missing or empty"
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
tryCreateSnpMatrixMissingConsensusRaiseSampleError()
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
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when some consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingConsensusRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    cfsan_snp_pipeline data lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Run cfsan_snp_pipeline snp_matrix -- fail because of missing some, but not all consensus fasta files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    rm "$tempDir/samples/sample1/consensus.fasta"
    rm "$tempDir/samples/sample4/consensus.fasta"
    rm "$tempDir/snpma.fasta"
    cfsan_snp_pipeline snp_matrix -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_matrix error handling behavior
    assertEquals "cfsan_snp_pipeline snp_matrix returned incorrect error code when some consensus.fasta missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingConsensusRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
tryCreateSnpReferenceSeqPermissionTrap()
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
    cfsan_snp_pipeline snp_reference -l "$tempDir/snplist.txt" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline snp_reference"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Error detected while running cfsan_snp_pipeline snp_reference"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/referenceSNP.fasta"
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqPermissionTrap 100
}

# Verify the cfsan_snp_pipeline snp_reference script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqPermissionTrap 100
}


# Verify the cfsan_snp_pipeline snp_reference script detects failure.
tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError()
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
    cfsan_snp_pipeline snp_reference -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/create_snp_reference_seq.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the snplist file"
    assertFileContains "$logDir/create_snp_reference_seq.log" "cannot create the snp reference sequence without the snplist file"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline snp_reference script detects failure.
tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError()
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
    cfsan_snp_pipeline snp_reference -l "$tempDir/snplist" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline snp_reference error handling behavior
    assertEquals "cfsan_snp_pipeline snp_reference returned incorrect error code when reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "cfsan_snp_pipeline snp_reference failed."
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/create_snp_reference_seq.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileContains "$logDir/create_snp_reference_seq.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline snp_reference script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
tryCollectSampleMetricsMissingSampleDirRaiseSampleError()
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
    cfsan_snp_pipeline collect_metrics -c "$tempDir/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metricsn error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the sample directory was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileContains "$logDir/collectSampleMetrics.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 98
}

# Verify the cfsan_snp_pipeline collect_metrics script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 100
}


# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
tryCollectSampleMetricsMissingReferenceRaiseGlobalError()
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
    cfsan_snp_pipeline collect_metrics -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/collectSampleMetrics.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
tryCollectSampleMetricsMissingInputFiles()
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
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."

    assertFileContains "$logDir/collectSampleMetrics.log" "SAM file reads.sam was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Deduped BAM file reads.sorted.deduped.bam was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "BAM file reads.sorted.bam was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Pileup file reads.all.pileup was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "VCF file var.flt.vcf was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "VCF file var.flt_preserved.vcf was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus VCF file consensus.vcf was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus VCF file consensus_preserved.vcf was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus fasta file consensus.fasta was not found"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus fasta file consensus_preserved.fasta was not found"

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
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Cannot calculate"
    assertFileContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectSampleMetricsMissingInputFilesStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsMissingInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectSampleMetricsMissingInputFilesNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsMissingInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles missing input files
testCollectSampleMetricsMissingInputFilesStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsMissingInputFiles 0
}


# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
tryCollectSampleMetricsEmptyInputFiles()
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
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."

    assertFileContains "$logDir/collectSampleMetrics.log" "SAM file reads.sam is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Deduped BAM file reads.sorted.deduped.bam is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "BAM file reads.sorted.bam is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Pileup file reads.all.pileup is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "VCF file var.flt.vcf is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "VCF file var.flt_preserved.vcf is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus VCF file consensus.vcf is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus VCF file consensus_preserved.vcf is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus fasta file consensus.fasta is empty"
    assertFileContains "$logDir/collectSampleMetrics.log" "Consensus fasta file consensus_preserved.fasta is empty"

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
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Cannot calculate"
    assertFileContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectSampleMetricsEmptyInputFilesStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsEmptyInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectSampleMetricsEmptyInputFilesNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsEmptyInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles empty input files
testCollectSampleMetricsEmptyInputFilesStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsEmptyInputFiles 0
}


# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
tryCollectSampleMetricsCorruptInputFiles()
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
    cfsan_snp_pipeline collect_metrics -o "$tempDir/samples/sample1/metrics" -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline collect_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline collect_metrics returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Error detected while running cfsan_snp_pipeline collect_metrics."
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
    assertFileContains "$logDir/collectSampleMetrics.log" "Cannot calculate number of reads and %mapped"
    assertFileContains "$tempDir/samples/sample1/metrics" "Cannot calculate mean insert size"
    assertFileContains "$logDir/collectSampleMetrics.log" "Cannot calculate mean insert size"
    assertFileContains "$tempDir/samples/sample1/metrics" "Cannot calculate mean pileup depth"
    assertFileContains "$logDir/collectSampleMetrics.log" "Cannot calculate mean pileup depth"
    assertFileContains "$logDir/collectSampleMetrics.log" "cfsan_snp_pipeline collect_metrics finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectSampleMetricsCorruptInputFilesStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsCorruptInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectSampleMetricsCorruptInputFilesNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsCorruptInputFiles 0
}

# Verify the cfsan_snp_pipeline collect_metrics script handles corrupt input files
testCollectSampleMetricsCorruptInputFilesStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsCorruptInputFiles 0
}


# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
tryCombineSampleMetricsMissingSampleDirRaiseGlobalError()
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
    cfsan_snp_pipeline combine_metrics "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics failed."
    assertFileNotContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/combineSampleMetrics.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline combine_metrics script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning()
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
    cfsan_snp_pipeline combine_metrics -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline combine_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when the sample metrics file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics warning"
    assertFileNotContains "$tempDir/error.log" "cfsan_snp_pipeline combine_metrics failed"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics failed"
    assertFileContains "$tempDir/metrics.tsv" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$tempDir/error.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$logDir/combineSampleMetrics.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$tempDir/metrics.tsv" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$tempDir/error.log" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$logDir/combineSampleMetrics.log" "Sample metrics file $tempDir/samples/sample4/metrics is empty"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the cfsan_snp_pipeline combine_metrics script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}


# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
tryCombineSampleMetricsPermissionTrap()
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
    cfsan_snp_pipeline combine_metrics -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline combine_metrics error handling behavior
    assertEquals "cfsan_snp_pipeline combine_metrics returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline combine_metrics"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Error detected while running cfsan_snp_pipeline combine_metrics"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/metrics.tsv"
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsPermissionTrap 100
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsPermissionTrap 100
}

# Verify the cfsan_snp_pipeline combine_metrics script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsPermissionTrap 100
}


# Verify the cfsan_snp_pipeline distance script detects missing input file
tryCalculateSnpDistancesMissingInputRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run cfsan_snp_pipeline distance with missing snpma.fasta
    cfsan_snp_pipeline distance -p pp -m mm "$tempDir/snpma.fasta" &> "$logDir/calcSnpDistances.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when input snp matrix file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline distance failed."
    assertFileNotContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance failed"
    assertFileContains "$tempDir/error.log" "Error: cannot calculate sequence distances without the snp matrix file"
    assertFileContains "$logDir/calcSnpDistances.log" "Error: cannot calculate sequence distances without the snp matrix file"
    assertFileContains "$tempDir/error.log" "SNP matrix file $tempDir/snpma.fasta does not exist"
    assertFileContains "$logDir/calcSnpDistances.log" "SNP matrix file $tempDir/snpma.fasta does not exist"
    assertFileNotContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/calcSnpDistances.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testCalculateSnpDistancesMissingInputRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCalculateSnpDistancesMissingInputRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testCalculateSnpDistancesMissingInputRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCalculateSnpDistancesMissingInputRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects missing input file
testCalculateSnpDistancesMissingInputRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCalculateSnpDistancesMissingInputRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline distance script detects no output file options
tryCalculateSnpDistancesMissingOutputOptionsRaiseGlobalError()
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
    cfsan_snp_pipeline distance "$tempDir/snpma.fasta" &> "$logDir/calcSnpDistances.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when both output options were missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "cfsan_snp_pipeline distance failed."
    assertFileNotContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance failed"
    assertFileContains "$tempDir/error.log" "Error: no output file specified"
    assertFileContains "$logDir/calcSnpDistances.log" "Error: no output file specified"
    assertFileNotContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/calcSnpDistances.log" "Use the -f option to force a rebuild"
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testCalculateSnpDistancesMissingOutputOptionsRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCalculateSnpDistancesMissingOutputOptionsRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testCalculateSnpDistancesMissingOutputOptionsRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCalculateSnpDistancesMissingOutputOptionsRaiseGlobalError 100
}

# Verify the cfsan_snp_pipeline distance script detects no output file options
testCalculateSnpDistancesMissingOutputOptionsRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCalculateSnpDistancesMissingOutputOptionsRaiseGlobalError 100
}


# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
tryCalculateSnpDistancesPermissionTrap()
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
    cfsan_snp_pipeline distance -p "$tempDir/pairwise" "$tempDir/snpma.fasta" &> "$logDir/calcSnpDistances.log"
    errorCode=$?

    # Verify cfsan_snp_pipeline distance error handling behavior
    assertEquals "cfsan_snp_pipeline distance returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline distance"
    assertFileNotContains "$logDir/calcSnpDistances.log" "Error detected while running cfsan_snp_pipeline distance"
    assertFileContains "$tempDir/error.log" "IOError|PermissionError"
    assertFileContains "$logDir/calcSnpDistances.log" "IOError|PermissionError"
    assertFileNotContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"
    assertFileNotContains "$logDir/calcSnpDistances.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/pairwise"
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testCalculateSnpDistancesPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCalculateSnpDistancesPermissionTrap 100
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testCalculateSnpDistancesPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCalculateSnpDistancesPermissionTrap 100
}

# Verify the cfsan_snp_pipeline distance script traps attempts to write to unwritable file
testCalculateSnpDistancesPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCalculateSnpDistancesPermissionTrap 100
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
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    assertFileContains "$tempDir/error.log" "Configuration file $tempDir/not-exist.conf does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Configuration file $tempDir/not-exist.conf does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    assertFileContains "$tempDir/error.log" "The file of samples directories, not-exist-file, does not exist"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "The file of samples directories, not-exist-file, does not exist"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineValidateSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

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
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    assertFileContains "$tempDir/error.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running cfsan_snp_pipeline index_ref."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2-build"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "cfsan_snp_pipeline index_ref finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"

    assertFileContains "$tempDir/error.log" "See also the log files in directory $tempDir/logs"
    assertFileContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample2/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample3/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample4/reads.sam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.deduped.bam"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapAlignSampleToReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapNoStopAllFail()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
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
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample2/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample3/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample4/reads.sam"

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
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
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
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam"
    assertFileNotContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample2/reads.sam"
    assertFileNotContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample3/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample4/reads.sam"
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
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
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
    verifyNonEmptyReadableFile "$logDir/prepReference.log"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-1"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-2"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-4"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-3"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-1"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-2"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-4"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-3"
    verifyNonEmptyReadableFile "$logDir/snpList.log"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-1"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-2"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-4"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-3"
    verifyNonEmptyReadableFile "$logDir/snpMatrix.log"
    verifyNonEmptyReadableFile "$logDir/snpReference.log"
    verifyNonEmptyReadableFile "$logDir/mergeVcf.log"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-1"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-2"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-4"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-3"
    verifyNonEmptyReadableFile "$logDir/combineSampleMetrics.log"
    verifyNonEmptyReadableFile "$logDir/calcSnpDistances.log"
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
    assertFileContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/alignSamples.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/alignSamples.log-2" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/alignSamples.log-3" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/alignSamples.log-4" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/prepSamples.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/prepSamples.log-2" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/prepSamples.log-3" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/prepSamples.log-4" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/snpList.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-3" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus.log-4" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectSampleMetrics.log-2" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectSampleMetrics.log-4" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/collectSampleMetrics.log-3" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/snpList_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-3" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/callConsensus_preserved.log-4" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/calcSnpDistances_preserved.log" "cfsan_snp_pipeline distance finished"

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
    assertFileContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/alignSamples.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/prepSamples.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/snpList.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/snpList_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/calcSnpDistances_preserved.log" "cfsan_snp_pipeline distance finished"

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
    touch -d '-10 day' $tempDir/reference/*.fasta
    touch -d  '-9 day' $tempDir/reference/*.bt2
    touch -d  '-8 day' $tempDir/samples/*/*.fastq
    touch -d  '-7 day' $tempDir/samples/*/reads.sam
    touch -d  '-6 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-4 day' $tempDir/samples/*/reads.sorted.deduped.bam
    touch -d  '-3 day' $tempDir/samples/*/reads.all.pileup
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
    assertFileContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/alignSamples.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/prepSamples.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/snpList.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/snpList_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/calcSnpDistances_preserved.log" "cfsan_snp_pipeline distance finished"
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
    touch -d  '-3 day' $tempDir/samples/*/reads.all.pileup
    touch $tempDir/samples/*/var.flt.vcf
    rm -rf "$tempDir/samples/sample1"
    sleep 1

    cfsan_snp_pipeline data configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

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
    assertFileContains "$logDir/prepSamples.log-2" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/snpList.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectSampleMetrics.log-2" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/snpList_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-2" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/calcSnpDistances_preserved.log" "cfsan_snp_pipeline distance finished"
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
    touch -d '-12 day' $tempDir/reference/*.fasta
    touch -d '-11 day' $tempDir/reference/*.bt2
    touch -d '-10 day' $tempDir/samples/*/*.fastq
    touch -d  '-9 day' $tempDir/samples/*/reads.sam
    touch -d  '-8 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-7 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-6 day' $tempDir/samples/*/reads.sorted.deduped.bam
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
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=23.22"

    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=20000"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberDupReads=110"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=94.55"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=286.84"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=23.22"
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
    verifyNonEmptyReadableFile "$logDir/prepReference.log"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-1"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-2"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-4"
    verifyNonEmptyReadableFile "$logDir/alignSamples.log-3"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-1"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-2"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-4"
    verifyNonEmptyReadableFile "$logDir/prepSamples.log-3"
    verifyNonEmptyReadableFile "$logDir/snpList.log"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-1"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-2"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-4"
    verifyNonEmptyReadableFile "$logDir/callConsensus.log-3"
    verifyNonEmptyReadableFile "$logDir/snpMatrix.log"
    verifyNonEmptyReadableFile "$logDir/snpReference.log"
    verifyNonEmptyReadableFile "$logDir/mergeVcf.log"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-1"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-2"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-4"
    verifyNonEmptyReadableFile "$logDir/collectSampleMetrics.log-3"

    assertFileContains "$logDir/prepReference.log" "lambda_virus.rev.1.bt2 is already freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepReference.log" "lambda_virus.fasta.fai is already freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/alignSamples.log-1" "sample1 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/alignSamples.log-2" "sample2 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/alignSamples.log-3" "sample4 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/alignSamples.log-4" "sample3 has already been aligned to lambda_virus.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/prepSamples.log-1" "Unsorted bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-2" "Unsorted bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-3" "Unsorted bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-4" "Unsorted bam file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-1" "Sorted bam file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-2" "Sorted bam file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-3" "Sorted bam file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-4" "Sorted bam file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-1" "Pileup file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-2" "Pileup file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-3" "Pileup file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-4" "Pileup file is already freshly created for sample3.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-1" "VCF file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-2" "VCF file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-3" "VCF file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-4" "VCF file is already freshly created for sample3.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpList.log" "snplist.txt has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/callConsensus.log-1" "sample1/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-2" "sample2/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-3" "sample4/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus.log-4" "sample3/consensus.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpMatrix.log" "/snpma.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpReference.log" "referenceSNP.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeVcf.log" "Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."

    assertFileNotContains "$logDir/collectSampleMetrics.log-1" "already freshly created"
    assertFileNotContains "$logDir/collectSampleMetrics.log-2" "already freshly created"
    assertFileNotContains "$logDir/collectSampleMetrics.log-3" "already freshly created"
    assertFileNotContains "$logDir/collectSampleMetrics.log-4" "already freshly created"

    assertFileContains "$logDir/calcSnpDistances.log" "have already been freshly built.  Use the -f option to force a rebuild"

    # =======
    assertFileContains "$logDir/filterAbnormalSNP.log" "All preserved and removed vcf files are already freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpList_preserved.log" "snplist_preserved.txt has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/callConsensus_preserved.log-1" "sample1/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-2" "sample2/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-3" "sample4/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/callConsensus_preserved.log-4" "sample3/consensus_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpMatrix_preserved.log" "/snpma_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/snpReference_preserved.log" "referenceSNP_preserved.fasta has already been freshly built.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/mergeVcf_preserved.log" "Multi-VCF file is already freshly created.  Use the -f option to force a rebuild."

    assertFileContains "$logDir/calcSnpDistances_preserved.log" "have already been freshly built.  Use the -f option to force a rebuild"

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
    echo 'CombineSampleMetrics_ExtraParams="-s"' >> "$tempDir/snppipeline.conf"
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
    sed -i s:SnpPipeline_MaxSnps=-1:SnpPipeline_MaxSnps=40: "$tempDir/snppipeline.conf"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

	# Verify each pipeline stage runs to completion
    logDir=$(echo $(ls -d $tempDir/logs*))
    assertFileContains "$logDir/prepReference.log" "cfsan_snp_pipeline index_ref finished"
    assertFileContains "$logDir/alignSamples.log-1" "cfsan_snp_pipeline map_reads finished"
    assertFileContains "$logDir/prepSamples.log-1" "cfsan_snp_pipeline call_sites finished"
    assertFileContains "$logDir/snpList.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "cfsan_snp_pipeline collect_metrics finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "cfsan_snp_pipeline combine_metrics finished"
    assertFileContains "$logDir/calcSnpDistances.log" "cfsan_snp_pipeline distance finished"

    assertFileContains "$logDir/filterAbnormalSNP.log" "cfsan_snp_pipeline filter_regions finished"
    assertFileContains "$logDir/snpList_preserved.log" "cfsan_snp_pipeline merge_sites finished"
    assertFileContains "$logDir/callConsensus_preserved.log-1" "cfsan_snp_pipeline call_consensus finished"
    assertFileContains "$logDir/mergeVcf_preserved.log" "cfsan_snp_pipeline merge_vcfs finished"
    assertFileContains "$logDir/snpMatrix_preserved.log" "cfsan_snp_pipeline snp_matrix finished"
    assertFileContains "$logDir/snpReference_preserved.log" "cfsan_snp_pipeline snp_reference finished"
    assertFileContains "$logDir/calcSnpDistances_preserved.log" "cfsan_snp_pipeline distance finished"

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
