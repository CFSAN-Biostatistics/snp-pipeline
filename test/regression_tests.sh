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
    grepResult=$(grep "$targetString" "$file")
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
    assertNotNull "copy_snppipeline_data.py is not on the path"     "$(which copy_snppipeline_data.py)"
    assertNotNull "run_snp_pipeline.sh is not on the path"          "$(which run_snp_pipeline.sh)"
    assertNotNull "prepReference.sh is not on the path"             "$(which prepReference.sh)"
    assertNotNull "alignSampleToReference.sh is not on the path"    "$(which alignSampleToReference.sh)"
    assertNotNull "prepSamples.sh is not on the path"               "$(which prepSamples.sh)"
    assertNotNull "call_consensus.py is not on the path"            "$(which call_consensus.py)"
    assertNotNull "create_snp_list.py is not on the path"           "$(which create_snp_list.py)"
    assertNotNull "create_snp_pileup.py is not on the path"         "$(which create_snp_pileup.py)"
    assertNotNull "create_snp_matrix.py is not on the path"         "$(which create_snp_matrix.py)"
    assertNotNull "create_snp_reference_seq.py is not on the path"  "$(which create_snp_reference_seq.py)"
    assertNotNull "collectSampleMetrics.sh is not on the path"      "$(which collectSampleMetrics.sh)"
    assertNotNull "combineSampleMetrics.sh is not on the path"      "$(which combineSampleMetrics.sh)"
}

# Verify copy_snppipeline_data.py emits lambda test data
testCopySnpPipelineLambdaData()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    copy_snppipeline_data.py lambdaVirusInputs $tempDir
    verifyNonEmptyReadableFile "$tempDir/reference/lambda_virus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/sample1_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/sample1_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/sample2_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample2/sample2_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/sample3_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample3/sample3_2.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/sample4_1.fastq"
    verifyNonEmptyReadableFile "$tempDir/samples/sample4/sample4_2.fastq"

    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
}

# Verify copy_snppipeline_data.py emits configuration file
testCopySnpPipelineConfigurationFile()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    copy_snppipeline_data.py configurationFile $tempDir
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
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "prepReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "alignSampleToReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "prepSamples.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_list.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "call_consensus.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_matrix.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "create_snp_reference_seq.py"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "collectSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "combineSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "$aligner"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "samtools"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "java"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "tabix"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bgzip"
    verifyWhetherCommandOnPathChecked "$tempDir/error.log" "bcftools"
    assertFileContains "$tempDir/error.log" "CLASSPATH is not configured with the path to VarScan"
    assertFileContains "$tempDir/error.log" "Check the SNP Pipeline installation instructions here: http://snp-pipeline.readthedocs.org/en/latest/installation.html"

    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "prepReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "alignSampleToReference.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "prepSamples.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_list.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "call_consensus.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_matrix.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "create_snp_reference_seq.py"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "collectSampleMetrics.sh"
    verifyWhetherCommandOnPathChecked "$tempDir/run_snp_pipeline.stderr.log" "combineSampleMetrics.sh"
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
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencyBowtie2RaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 bowtie2
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}

# Verify run_snp_pipeline checks for necessary scripts and tools
testRunSnpPipelineDependencySmaltRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    echo "SnpPipeline_Aligner=smalt" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineDependencyRaiseFatalError 1 smalt
}


# Verify the prepReference script detects a misconfigured environment variable
tryPrepReferenceEnvironmentRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables prepReference will use.
    # This simulates what run_snp_pipeline does before running prepReference.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately misconfigure the environment
    export SnpPipeline_Aligner=garbage

    # Run prepReference
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "prepReference.sh returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "prepReference.sh failed"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/prepReference.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEnvironmentRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceEnvironmentRaiseGlobalError 100
}

# Verify the prepReference script detects a misconfigured environment variable
tryPrepReferenceEmptyFastaFileRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables prepReference will use.
    # This simulates what run_snp_pipeline does before running prepReference.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately create empty fasta file
    rm "$tempDir/reference/lambda_virus.fasta"
    touch "$tempDir/reference/lambda_virus.fasta"

    # Run prepReference
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "prepReference.sh returned incorrect error code when the reference file was empty." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "prepReference.sh failed"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileContains "$logDir/prepReference.log" "Reference file $tempDir/reference/lambda_virus.fasta is empty"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}

# Verify the prepReference script detects a misconfigured environment variable
testPrepReferenceEmptyFastaFileRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceEmptyFastaFileRaiseGlobalError 100
}


# Verify the prepReference script detects bowtie error and emits the global error marker file.
tryPrepReferenceBowtieIndexTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables prepReference will use.
    # This simulates what run_snp_pipeline does before running prepReference.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"

    # Run prepReference
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "prepReference.sh / bowtie returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepReference.sh."
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running prepReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the prepReference script detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceBowtieIndexTrap 100
}

# Verify the prepReference script detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceBowtieIndexTrap 100
}

# Verify the prepReference script detects bowtie error and emits the global error marker file.
testPrepReferenceBowtieIndexTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceBowtieIndexTrap 100
}


# Verify the prepReference script detects smalt error and emits the global error marker file.
tryPrepReferenceSmaltIndexTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables prepReference will use.
    # This simulates what run_snp_pipeline does before running prepReference.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"
    sed -i 's/A/>\n/g' "$tempDir/reference/lambda_virus.fasta"

    # Run prepReference
    export SnpPipeline_Aligner=smalt
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "prepReference.sh / smalt returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepReference.sh."
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running prepReference.sh."
    assertFileContains "$tempDir/error.log" "smalt index"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the prepReference script detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the prepReference script detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the prepReference script detects smalt error and emits the global error marker file.
testPrepReferenceSmaltIndexTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceSmaltIndexTrap 100
}

# Verify the prepReference script detects samtools error and emits the global error marker file.
tryPrepReferenceSamtoolsFaidxTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables prepReference will use.
    # This simulates what run_snp_pipeline does before running prepReference.
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt the fasta file, changing line lengths
    sed -i 's/A/AA/g' "$tempDir/reference/lambda_virus.fasta"

    # Run prepReference
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Verify error handling behavior
    assertEquals "prepReference.sh / samtools returned incorrect error code when the input fasta was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepReference.sh."
    assertFileNotContains "$logDir/prepReference.log" "Error detected while running prepReference.sh."
    assertFileContains "$tempDir/error.log" "samtools faidx"
    assertFileNotContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileNotContains "$logDir/prepReference.log" "Use the -f option to force a rebuild"
}

# Verify the prepReference script detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepReferenceSamtoolsFaidxTrap 100
}

# Verify the prepReference script detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepReferenceSamtoolsFaidxTrap 100
}

# Verify the prepReference script detects samtools error and emits the global error marker file.
testPrepReferenceSamtoolsFaidxTrapUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepReferenceSamtoolsFaidxTrap 100
}



# Verify the alignSampleToReference script detects a misconfigured environment variable.
tryAlignSampleToReferenceEnvironmentRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately misconfigure the environment
    export SnpPipeline_Aligner=garbage

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference returned incorrect error code when the SnpPipeline_Aligner environment variable was misconfigured." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh failed"
    assertFileContains "$tempDir/error.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileContains "$logDir/alignSampleToReference.log" "Error: only bowtie2 and smalt aligners are supported."
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the alignSampleToReference script detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}

# Verify the alignSampleToReference script detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}

# Verify the alignSampleToReference script detects a misconfigured environment variable.
testAlignSampleToReferenceEnvironmentRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceEnvironmentRaiseGlobalError 100
}


# Verify the alignSampleToReference script detects an Missing reference file
tryAlignSampleToReferenceMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"
}

# Verify the alignSampleToReference script detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the alignSampleToReference script detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}

# Verify the alignSampleToReference script detects a missing reference file
testAlignSampleToReferenceMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingReferenceRaiseGlobalError 100
}


# Verify the alignSampleToReference script detects a missing sample file
tryAlignSampleToReferenceMissingSample1RaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Sample file $tempDir/samples/sample1/sample1_1.fastq does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 100
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 98
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample1 variable.
testAlignSampleToReferenceMissingSample1RaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingSample1RaiseSampleError 100
}


# Verify the alignSampleToReference script detects a missing sample file
tryAlignSampleToReferenceMissingSample2RaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try align with missing file
    rm "$tempDir/samples/sample1/sample1_2.fastq"

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference returned incorrect error code when the sample file 1 was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh failed"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh failed"
    assertFileContains "$tempDir/error.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileContains "$logDir/alignSampleToReference.log" "Sample file $tempDir/samples/sample1/sample1_2.fastq does not exist"
    assertFileNotContains "$logDir/alignSampleToReference.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 100
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 98
}

# Verify the alignSampleToReference script detects a misconfigured MissingSample2 variable.
testAlignSampleToReferenceMissingSample2RaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceMissingSample2RaiseSampleError 100
}


# Verify the alignSampleToReference script detects bowtie alignment error.
tryAlignSampleToReferenceBowtieAlignTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"
    export SnpPipeline_Aligner=bowtie2

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log-1"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/alignSampleToReference.log-1" "bowtie2"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/alignSampleToReference.log-3"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference with bowtie2 returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$logDir/alignSampleToReference.log-3" "bowtie2"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Use the -f option to force a rebuild"
}


# Verify the alignSampleToReference script detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceBowtieAlignTrap 100
}

# Verify the alignSampleToReference script detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceBowtieAlignTrap 98
}

# Verify the alignSampleToReference script detects bowtie alignment error.
testAlignSampleToReferenceBowtieAlignTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceBowtieAlignTrap 100
}


# Verify the alignSampleToReference script detects smalt alignment error.
tryAlignSampleToReferenceSmaltAlignTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"
    export SnpPipeline_Aligner=smalt

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"

    # Deliberately corrupt the FASTQ file
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"

    # Try to align a paired-end sample to the reference
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> "$logDir/alignSampleToReference.log-1"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/alignSampleToReference.log-1" "smalt map"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-1" "Use the -f option to force a rebuild"

    # Repeat the test with an unpaired fastq
    echo "Garbage" > "$tempDir/samples/sample3/sample3_1.fastq"
    rm "$tempDir/error.log" &> /dev/null
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample3/sample3_1.fastq" &> "$logDir/alignSampleToReference.log-3"
    errorCode=$?

    # Verify alignSampleToReference error handling behavior
    assertEquals "alignSampleToReference with smalt returned incorrect error code when the input fastq was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "smalt map"
    assertFileContains "$logDir/alignSampleToReference.log-3" "smalt map"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "alignSampleToReference.sh finished"
    assertFileNotContains "$logDir/alignSampleToReference.log-3" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}


# Verify the alignSampleToReference script detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryAlignSampleToReferenceSmaltAlignTrap 100
}

# Verify the alignSampleToReference script detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryAlignSampleToReferenceSmaltAlignTrap 98
}

# Verify the alignSampleToReference script detects smalt alignment error.
testAlignSampleToReferenceSmaltAlignTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryAlignSampleToReferenceSmaltAlignTrap 100
}


# Verify the PrepSamples script detects a Missing reference file
tryPrepSamplesMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try prepSamples with missing file
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to align a paired-end sample to the reference
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta" "xxxx" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify PrepSamples error handling behavior
    assertEquals "PrepSamples returned incorrect error code when the reference file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "prepSamples.sh failed"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/prepSamples.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the PrepSamples script detects a missing Reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}

# Verify the PrepSamples script detects a missing reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}

# Verify the PrepSamples script detects a missing reference file
testPrepSamplesMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesMissingReferenceRaiseGlobalError 100
}


# Verify the PrepSamples script detects a missing sample sam file
tryPrepSamplesMissingSamFileRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately try prepSamples with missing file
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify PrepSamples error handling behavior
    assertEquals "PrepSamples returned incorrect error code when the reads.sam file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "prepSamples.sh failed"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh failed"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam does not exist"
    assertFileContains "$logDir/prepSamples.log" "Sample SAM file $tempDir/samples/sample1/reads.sam does not exist"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"

    # Restore the normal aligner
    export SnpPipeline_Aligner=bowtie2
}

# Verify the PrepSamples script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesMissingSamFileRaiseSampleError 100
}

# Verify the PrepSamples script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesMissingSamFileRaiseSampleError 98
}

# Verify the PrepSamples script detects a misconfigured MissingSamFile variable.
testPrepSamplesMissingSamFileRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesMissingSamFileRaiseSampleError 100
}


# Verify the prepSamples script detects Samtools view failure.
tryPrepSamplesSamtoolsViewTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run prepSamples normally
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the SAM file and re-run prepSamples.sh
    echo "Garbage" > "$tempDir/samples/sample1/reads.sam"
    rm "$tempDir/error.log" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify prepSamples error handling behavior
    assertEquals "prepSamples.sh returned incorrect error code when the input SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepSamples.sh."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running prepSamples.sh."
    assertFileContains "$tempDir/error.log" "samtools view"
    assertFileContains "$logDir/prepSamples.log" "samtools view"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Use the -f option to force a rebuild"
}

# Verify the prepSamples script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsViewTrap 100
}

# Verify the prepSamples script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsViewTrap 98
}

# Verify the prepSamples script detects Samtools view failure.
testPrepSamplesSamtoolsViewTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsViewTrap 100
}


# Verify the prepSamples script detects Samtools sort failure.
tryPrepSamplesSamtoolsSortTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.  
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run prepSamples normally
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the reads.unsorted.bam file and re-run prepSamples.sh
    echo "Garbage" > "$tempDir/samples/sample1/reads.unsorted.bam"
    rm "$tempDir/error.log" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify prepSamples error handling behavior
    assertEquals "prepSamples.sh returned incorrect error code when reads.unsorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepSamples.sh."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running prepSamples.sh."
    assertFileContains "$tempDir/error.log" "samtools sort"
    assertFileContains "$logDir/prepSamples.log" "samtools sort"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Sorted bam file is already freshly created"
}

# Verify the prepSamples script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsSortTrap 100
}

# Verify the prepSamples script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsSortTrap 98
}

# Verify the prepSamples script detects Samtools sort failure.
testPrepSamplesSamtoolsSortTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsSortTrap 100
}


# Verify the prepSamples script detects Samtools mpileup failure.
tryPrepSamplesSamtoolsMpileupTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run prepSamples normally
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Deliberately corrupt the reads.sorted.bam file and re-run prepSamples.sh
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    echo "Garbage" > "$tempDir/samples/sample1/reads.sorted.bam"
    rm "$tempDir/error.log" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?

    # Verify prepSamples error handling behavior
    assertEquals "prepSamples.sh returned incorrect error code when reads.sorted.bam was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running prepSamples.sh."
    assertFileNotContains "$logDir/prepSamples.log" "Error detected while running prepSamples.sh."
    assertFileContains "$tempDir/error.log" "samtools mpileup"
    assertFileContains "$logDir/prepSamples.log" "samtools mpileup"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Pileup file is already freshly created"
}

# Verify the prepSamples script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesSamtoolsMpileupTrap 100
}

# Verify the prepSamples script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesSamtoolsMpileupTrap 98
}

# Verify the prepSamples script detects Samtools mpileup failure.
testPrepSamplesSamtoolsMpileupTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesSamtoolsMpileupTrap 100
}


# Verify the prepSamples script detects Varscan failure.
tryPrepSamplesVarscanRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run prepSamples normally
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Configure VarScan with invalid parameter settings and re-run prepSamples.sh
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.all.pileup"
    export VarscanMpileup2snp_ExtraParams="--min-coverage -1 --min-reads 99999999 --min-avg_qual -100 --min-var-freq 2 --output-vcf 2 --invalid-parameter"
    rm "$tempDir/error.log" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?
    export VarscanMpileup2snp_ExtraParams=""

    # Verify prepSamples error handling behavior
    assertEquals "prepSamples.sh returned incorrect error code when varscan failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/prepSamples.log" "var.flt.vcf is empty"
    assertFileContains "$logDir/prepSamples.log" "VarScan mpileup2snp"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Vcf file is already freshly created"
}

# Verify the prepSamples script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesVarscanRaiseSampleError 100
}

# Verify the prepSamples script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesVarscanRaiseSampleError 98
}

# Verify the prepSamples script detects Varscan failure.
testPrepSamplesVarscanRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesVarscanRaiseSampleError 100
}


# Verify the prepSamples script detects unset java classpath needed to run Varscan.
tryPrepSamplesVarscanClasspathRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    export Bowtie2Align_ExtraParams="--reorder -X 1000"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null

    # First run prepSamples normally
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/var.flt.vcf"

    # Configure VarScan with invalid parameter settings and re-run prepSamples.sh
    touch "$tempDir/samples/sample1/reads.unsorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.sorted.bam"
    sleep 1
    touch "$tempDir/samples/sample1/reads.all.pileup"
    saveClassPath="$CLASSPATH"
    unset CLASSPATH
    rm "$tempDir/error.log" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"
    errorCode=$?
    export CLASSPATH="$saveClassPath"

    # Verify prepSamples error handling behavior
    assertEquals "prepSamples.sh returned incorrect error code when CLASSPATH unset failed." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "prepSamples.sh failed"
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh failed"
    assertFileContains "$tempDir/error.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileContains "$logDir/prepSamples.log" "Error: cannot execute VarScan. Define the path to VarScan in the CLASSPATH environment variable."
    assertFileNotContains "$logDir/prepSamples.log" "prepSamples.sh finished"
    assertFileNotContains "$logDir/prepSamples.log" "Vcf file is already freshly created"
}

# Verify the prepSamples script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the prepSamples script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}

# Verify the prepSamples script detects unset java classpath needed to run Varscan.
testPrepSamplesVarscanClasspathRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryPrepSamplesVarscanClasspathRaiseGlobalError 100 # this is a global error because all samples will fail
}


# Verify the create_snp_list.py script traps attempts to write to unwritable file
tryCreateSnpListPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    create_snp_list.py -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify create_snp_list.py error handling behavior
    assertEquals "create_snp_list.py returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running create_snp_list.py"
    assertFileNotContains "$logDir/create_snp_list.log" "Error detected while running create_snp_list.py"
    assertFileContains "$tempDir/error.log" "IOError"
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snplist.txt"
}

# Verify the create_snp_list.py script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListPermissionTrap 100
}

# Verify the create_snp_list.py script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListPermissionTrap 100
}

# Verify the create_snp_list.py script traps attempts to write to unwritable file
testCreateSnpListPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListPermissionTrap 100
}



# Verify the create_snp_list script detects failure.
tryCreateSnpListMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_list.py with missing sampleDirectories.txt
    create_snp_list.py -o "$tempDir/snplist.txt"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify create_snp_list error handling behavior
    assertEquals "create_snp_list.py returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_list.py failed."
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/create_snp_list.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingSampleDirRaiseGlobalError 100
}


# Verify the create_snp_list script detects failure.
tryCreateSnpListMissingVcfRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_list.py -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    create_snp_list.py -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify create_snp_list error handling behavior
    assertEquals "create_snp_list.py returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_list.py failed."
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py failed."
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
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingVcfRaiseGlobalError 100
}


# Verify the create_snp_list script detects failure.
tryCreateSnpListMissingVcfRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run create_snp_list.py -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    create_snp_list.py -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify create_snp_list error handling behavior
    assertEquals "create_snp_list.py returned incorrect error code when var.flt.vcf was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_list.py failed."
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py failed."
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpListMissingVcfRaiseSampleError 100
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run create_snp_list.py -- fail because of missing some, but not all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    create_snp_list.py -o "$tempDir/snplist.txt"  "$tempDir/sampleDirList.txt" &> "$logDir/create_snp_list.log"
    errorCode=$?

    # Verify create_snp_list error handling behavior
    assertEquals "create_snp_list.py returned incorrect error code when var.flt.vcf was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_list.py"
    assertFileNotContains "$tempDir/error.log" "create_snp_list.py failed."
    assertFileNotContains "$logDir/create_snp_list.log" "create_snp_list.py failed."
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample2/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample3/var.flt.vcf does not exist"
    assertFileContains "$logDir/create_snp_list.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "Error: 3 VCF files were missing or empty"
    assertFileContains "$logDir/create_snp_list.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/create_snp_list.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_list script detects failure.
testCreateSnpListMissingVcfRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpListMissingVcfRaiseSampleError 100
}


# Verify the call_consensus script traps on corrupt snplist file
tryCallConsensusCorruptSnplistTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run call_consensus.py -- fail because of corrupt snplist
    echo "Corrupt snplist content" > "$tempDir/snplist.txt"
    echo "Dummy pileup content" > "$tempDir/samples/sample1/reads.all.pileup"
    call_consensus.py -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "call_consensus.py returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running call_consensus.py"
    assertFileNotContains "$logDir/call_consensus.log" "Error detected while running call_consensus.py"
    assertFileContains "$tempDir/error.log" "function read_snp_position_list at line"
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus script traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusCorruptSnplistTrap 100
}

# Verify the call_consensus script traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusCorruptSnplistTrap 98
}

# Verify the call_consensus script traps on corrupt snplist file
testCallConsensusCorruptSnplistTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCallConsensusCorruptSnplistTrap 100
}




# Verify the call_consensus script detects failure.
tryCallConsensusMissingSnpListRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run call_consensus.py -- fail because of missing snplist
    call_consensus.py -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify call_consensus error handling behavior
    assertEquals "call_consensus.py returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "call_consensus.py failed."
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the snplist file"
    assertFileContains "$logDir/call_consensus.log" "cannot call consensus without the snplist file"
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/call_consensus.log" "Snplist file snplist.txt does not exist"
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCallConsensusMissingSnpListRaiseGlobalError 100
}


# Verify the call_consensus script does not fail merely because the snplist file is empty.
tryCallConsensusEmptySnpList()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run call_consensus.py with empty snplist
    touch "$tempDir/snplist.txt"
    call_consensus.py -l "$tempDir/snplist.txt" -o "$tempDir/samples/sample1/consensus.fasta"  --vcfFileName "$tempDir/samples/sample1/consensus.vcf" "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify call_consensus error handling behavior
    assertEquals "call_consensus.py returned incorrect error code when snplist was empty." $expectErrorCode $errorCode
    verifyNonExistingFile "$tempDir/error.log"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonEmptyReadableFile "$tempDir/samples/sample1/consensus.vcf"
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py failed."
    assertFileNotContains "$logDir/call_consensus.log" "cannot call consensus without the snplist file"
    assertFileNotContains "$logDir/call_consensus.log" "Snplist file $tempDir/snplist.txt is empty"
    assertFileContains "$logDir/call_consensus.log" "call_consensus.py finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus script does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus script does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusEmptySnpList 0
}

# Verify the call_consensus script does not fail merely because the snplist file is empty.
testCallConsensusEmptySnpListStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCallConsensusEmptySnpList 0
}


# Verify the call_consensus script detects failure.
tryCallConsensusMissingPileupRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    echo "fake snplist" > $tempDir/snplist.txt

    # Run call_consensus.py -- fail because of missing pileup
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    call_consensus.py -l "$tempDir/snplist.txt" -o "$tempDir/consensus.fasta"  "$tempDir/samples/sample1/reads.all.pileup" &> "$logDir/call_consensus.log"
    errorCode=$?

    # Verify call_consensus error handling behavior
    assertEquals "call_consensus.py returned incorrect error code when pileup file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "call_consensus.py failed."
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py failed."
    assertFileContains "$tempDir/error.log" "cannot call consensus without the pileup file"
    assertFileContains "$logDir/call_consensus.log" "cannot call consensus without the pileup file"
    assertFileContains "$tempDir/error.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileContains "$logDir/call_consensus.log" "Pileup file $tempDir/samples/sample1/reads.all.pileup does not exist"
    assertFileNotContains "$logDir/call_consensus.log" "call_consensus.py finished"
    assertFileNotContains "$logDir/call_consensus.log" "Use the -f option to force a rebuild"
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingPileupRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCallConsensusMissingPileupRaiseSampleError 100
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingPileupRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCallConsensusMissingPileupRaiseSampleError 98
}

# Verify the call_consensus script detects failure.
testCallConsensusMissingPileupRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCallConsensusMissingPileupRaiseSampleError 100
}


# Verify the mergeVcf script detects failure.
tryMergeVcfCorruptVcfTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run mergeVcf.sh with corrupt consensus.vcf
    sleep 1
    sed -i 's/1/@@@/g' "$tempDir/samples/sample1/consensus.vcf"
    mergeVcf.sh -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify mergeVcf error handling behavior
    assertEquals "mergeVcf.sh returned incorrect error code when consensus.vcf was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running mergeVcf.sh."
    assertFileNotContains "$logDir/mergeVcf.log" "Error detected while running mergeVcf.sh."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the mergeVcf script detects failure.
testMergeVcfCorruptVcfTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfCorruptVcfTrap 100
}

# Verify the mergeVcf script detects failure.
testMergeVcfCorruptVcfTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfCorruptVcfTrap 100
}

# Verify the mergeVcf script detects failure.
testMergeVcfCorruptVcfTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfCorruptVcfTrap 100
}


# Verify the mergeVcf script detects failure.
tryMergeVcfMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run mergeVcf.sh with missing sampleDirectories.txt
    mergeVcf.sh -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify mergeVcf error handling behavior
    assertEquals "mergeVcf.sh returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "mergeVcf.sh failed."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh failed."
    assertFileContains "$tempDir/error.log" "Sample directories file $tempDir/sampleDirectories.txt does not exist."
    assertFileContains "$logDir/mergeVcf.log" "Sample directories file $tempDir/sampleDirectories.txt does not exist."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfMissingSampleDirRaiseGlobalError 100
}


# Verify the mergeVcf script detects failure.
tryMergeVcfMissingVcfRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run mergeVcf.sh with missing vcf files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    mergeVcf.sh -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify mergeVcf error handling behavior
    assertEquals "mergeVcf.sh returned incorrect error code when vcf file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "mergeVcf.sh failed."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingVcfRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryMergeVcfMissingVcfRaiseSampleError 100
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingVcfRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run mergeVcf.sh with missing vcf files
    rm "$tempDir/samples/sample1/consensus.vcf"
    rm "$tempDir/snpma.vcf"
    mergeVcf.sh -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify mergeVcf keeps running when only one vcf file is missing
    assertEquals "mergeVcf.sh returned incorrect error code when vcf file was missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "mergeVcf.sh"
    assertFileNotContains "$tempDir/error.log" "mergeVcf.sh failed."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh failed."
    assertFileContains "$tempDir/error.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "Sample vcf file $tempDir/samples/sample1/consensus.vcf does not exist."
    assertFileContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
}

# Verify the mergeVcf script detects failure.
testMergeVcfMissingVcfRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryMergeVcfMissingVcfRaiseSampleError 100
}


# Verify the mergeVcf script detects failure.
tryMergeVcfZeroGoodSamplesRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run mergeVcf.sh with no consensus vcf files
    echo   "$tempDir/samples/sample1" > "$tempDir/sampleDirectories.txt"
    echo   "$tempDir/samples/sample2" >> "$tempDir/sampleDirectories.txt"
    mergeVcf.sh -o "$tempDir/snpma.vcf"  "$tempDir/sampleDirectories.txt" &> "$logDir/mergeVcf.log"
    errorCode=$?

    # Verify mergeVcf error handling behavior
    assertEquals "mergeVcf.sh returned incorrect error code when no good VCF files." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "mergeVcf.sh failed."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh failed."
    assertFileContains "$tempDir/error.log" "There are no vcf files to merge."
    assertFileContains "$logDir/mergeVcf.log" "There are no vcf files to merge."
    assertFileNotContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileNotContains "$logDir/mergeVcf.log" "Use the -f option to force a rebuild"
}

# Verify the mergeVcf script detects failure.
#testMergeVcfZeroGoodSamplesRaiseGlobalErrorStop()
#{
#    # Nothing to test, SnpPipeline_StopOnSampleError must be false to test this code path
#}

# Verify the mergeVcf script detects failure.
testMergeVcfZeroGoodSamplesRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryMergeVcfZeroGoodSamplesRaiseGlobalError 100
}

# Verify the mergeVcf script detects failure.
#testMergeVcfZeroGoodSamplesRaiseGlobalErrorStopUnset()
#{
#    # Nothing to test, SnpPipeline_StopOnSampleError must be false to test this code path
#}


# Verify the create_snp_matrix.py script traps attempts to write to unwritable file
tryCreateSnpMatrixPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    create_snp_matrix.py -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify create_snp_list.py error handling behavior
    assertEquals "create_snp_matrix.py returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running create_snp_matrix.py"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Error detected while running create_snp_matrix.py"
    assertFileContains "$tempDir/error.log" "IOError"
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/snpma.fasta"
}

# Verify the create_snp_matrix.py script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixPermissionTrap 100
}

# Verify the create_snp_matrix.py script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixPermissionTrap 100
}

# Verify the create_snp_matrix.py script traps attempts to write to unwritable file
testCreateSnpMatrixPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixPermissionTrap 100
}


# Verify the create_snp_matrix script detects failure.
tryCreateSnpMatrixMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_matrix.py with missing sampleDirectories.txt
    create_snp_matrix.py -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify create_snp_matrix error handling behavior
    assertEquals "create_snp_matrix.py returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_matrix.py failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py failed."
    assertFileContains "$tempDir/error.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "File of sample directories $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingSampleDirRaiseGlobalError 100
}


# Verify the create_snp_matrix script detects failure.
tryCreateSnpMatrixMissingConsensusRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run prep work
    #prepReference.sh "$tempDir/reference/lambda_virus.fasta" &> "$logDir/prepReference.log"
    #alignSampleToReference.sh "$tempDir/reference/lambda_virus.fasta" "$tempDir/samples/sample1/sample1_1.fastq" "$tempDir/samples/sample1/sample1_2.fastq" &> /dev/null
    #prepSamples.sh "$tempDir/reference/lambda_virus.fasta"  "$tempDir/samples/sample1" &> "$logDir/prepSamples.log"

    # Run create_snp_matrix.py -- fail because of missing all VCF files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirList.txt"
    create_snp_matrix.py -o "$tempDir/snpmap.fasta"  "$tempDir/sampleDirList.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify create_snp_matrix error handling behavior
    assertEquals "create_snp_matrix.py returned incorrect error code when all consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_matrix.py failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py failed."
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
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingConsensusRaiseGlobalError 100
}


# Verify the create_snp_matrix script detects failure.
tryCreateSnpMatrixMissingConsensusRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Run create_snp_matrix.py -- fail because of missing some, but not all consensus fasta files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    rm "$tempDir/samples/sample1/consensus.fasta"
    rm "$tempDir/samples/sample4/consensus.fasta"
    create_snp_matrix.py -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify create_snp_matrix error handling behavior
    assertEquals "create_snp_matrix.py returned incorrect error code when some consensus.fasta missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_matrix.py failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpMatrixMissingConsensusRaiseSampleError 100
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export SnpPipeline_StopOnSampleError=false
    export errorOutputFile="$tempDir/error.log"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Run create_snp_matrix.py -- fail because of missing some, but not all consensus fasta files
    printf "%s\n" $tempDir/samples/* >  "$tempDir/sampleDirectories.txt"
    rm "$tempDir/samples/sample1/consensus.fasta"
    rm "$tempDir/samples/sample4/consensus.fasta"
    rm "$tempDir/snpma.fasta"
    create_snp_matrix.py -o "$tempDir/snpma.fasta"  "$tempDir/sampleDirectories.txt" &> "$logDir/create_snp_matrix.log"
    errorCode=$?

    # Verify create_snp_matrix error handling behavior
    assertEquals "create_snp_matrix.py returned incorrect error code when some consensus.fasta missing." 0 $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_matrix.py"
    assertFileNotContains "$tempDir/error.log" "create_snp_matrix.py failed."
    assertFileNotContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py failed."
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample1/consensus.fasta does not exist"
    assertFileContains "$logDir/create_snp_matrix.log" "Consensus fasta file $tempDir/samples/sample4/consensus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "Error: 2 consensus fasta files were missing or empty"
    assertFileContains "$logDir/create_snp_matrix.log" "create_snp_matrix.py finished"
    assertFileNotContains "$logDir/create_snp_matrix.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_matrix script detects failure.
testCreateSnpMatrixMissingConsensusRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpMatrixMissingConsensusRaiseSampleError 100
}


# Verify the create_snp_reference_seq.py script traps attempts to write to unwritable file
tryCreateSnpReferenceSeqPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    create_snp_reference_seq.py -l "$tempDir/snplist.txt" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify create_snp_list.py error handling behavior
    assertEquals "create_snp_reference_seq.py returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running create_snp_reference_seq.py"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Error detected while running create_snp_reference_seq.py"
    assertFileContains "$tempDir/error.log" "IOError"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "create_snp_reference_seq.py finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/referenceSNP.fasta"
}

# Verify the create_snp_reference_seq.py script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqPermissionTrap 100
}

# Verify the create_snp_reference_seq.py script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqPermissionTrap 100
}

# Verify the create_snp_reference_seq.py script traps attempts to write to unwritable file
testCreateSnpReferenceSeqPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqPermissionTrap 100
}


# Verify the create_snp_reference_seq.py script detects failure.
tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_reference_seq.py -- fail because of missing snplist
    create_snp_reference_seq.py -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify create_snp_reference_seq.py error handling behavior
    assertEquals "create_snp_reference_seq.py returned incorrect error code when snplist was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_reference_seq.py failed."
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "create_snp_reference_seq.py failed."
    assertFileContains "$tempDir/error.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$logDir/create_snp_reference_seq.log" "Snplist file snplist.txt does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the snplist file"
    assertFileContains "$logDir/create_snp_reference_seq.log" "cannot create the snp reference sequence without the snplist file"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "create_snp_reference_seq.py finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingSnpListRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqMissingSnpListRaiseGlobalError 100
}


# Verify the create_snp_reference_seq.py script detects failure.
tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_reference_seq.py -- fail because of missing reference
    echo "Dummy snplist content" > "$tempDir/snplist"
    rm "$tempDir/reference/lambda_virus.fasta"
    create_snp_reference_seq.py -l "$tempDir/snplist" -o "$tempDir/referenceSNP.fasta"  "$tempDir/reference/lambda_virus.fasta" &> "$logDir/create_snp_reference_seq.log"
    errorCode=$?

    # Verify create_snp_reference_seq.py error handling behavior
    assertEquals "create_snp_reference_seq.py returned incorrect error code when reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "create_snp_reference_seq.py failed."
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "create_snp_reference_seq.py failed."
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/create_snp_reference_seq.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$tempDir/error.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileContains "$logDir/create_snp_reference_seq.log" "cannot create the snp reference sequence without the reference fasta file"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "create_snp_reference_seq.py finished"
    assertFileNotContains "$logDir/create_snp_reference_seq.log" "Use the -f option to force a rebuild"
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}

# Verify the create_snp_reference_seq.py script detects failure.
testCreateSnpReferenceSeqMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCreateSnpReferenceSeqMissingReferenceRaiseGlobalError 100
}


# Verify the collectSampleMetrics.sh script detects a missing sample directory
tryCollectSampleMetricsMissingSampleDirRaiseSampleError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately remove directory
    rm -rf "$tempDir/samples/sample1"

    # Try to collect metrics
    echo "Dummy consensus.fasta content" > "$tempDir/consensus.fasta"
    collectSampleMetrics.sh -c "$tempDir/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify collectSampleMetrics.sh error handling behavior
    assertEquals "collectSampleMetrics.sh returned incorrect error code when the sample directory was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "collectSampleMetrics.sh failed"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "collectSampleMetrics.sh failed"
    assertFileContains "$tempDir/error.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileContains "$logDir/collectSampleMetrics.log" "Sample directory $tempDir/samples/sample1 does not exist"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "collectSampleMetrics.sh finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the collectSampleMetrics.sh script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 100
}

# Verify the collectSampleMetrics.sh script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 98
}

# Verify the collectSampleMetrics.sh script detects a missing sample directory
testCollectSampleMetricsMissingSampleDirRaiseSampleErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsMissingSampleDirRaiseSampleError 100
}


# Verify the collectSampleMetrics.sh script detects missing reference
tryCollectSampleMetricsMissingReferenceRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately remove reference
    rm "$tempDir/reference/lambda_virus.fasta"

    # Try to collect metrics
    echo "Dummy consensus.fasta content" > "$tempDir/samples/sample1/consensus.fasta"
    collectSampleMetrics.sh -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify collectSampleMetrics.sh error handling behavior
    assertEquals "collectSampleMetrics.sh returned incorrect error code when the reference fasta file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "collectSampleMetrics.sh failed"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "collectSampleMetrics.sh failed"
    assertFileContains "$tempDir/error.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileContains "$logDir/collectSampleMetrics.log" "Reference file $tempDir/reference/lambda_virus.fasta does not exist"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "collectSampleMetrics.sh finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the collectSampleMetrics.sh script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the collectSampleMetrics.sh script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}

# Verify the collectSampleMetrics.sh script detects missing reference
testCollectSampleMetricsMissingReferenceRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsMissingReferenceRaiseGlobalError 100
}


# Verify the collectSampleMetrics.sh script traps errors
tryCollectSampleMetricsSamtoolsViewTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Deliberately corrupt files
    echo "Garbage" > "$tempDir/samples/sample1/reads.sam"

    # Try to collect metrics
    echo "Dummy consensus.fasta content" > "$tempDir/samples/sample1/consensus.fasta"
    collectSampleMetrics.sh -c "$tempDir/samples/sample1/consensus.fasta" "$tempDir/samples/sample1" "$tempDir/reference/lambda_virus.fasta" &> "$logDir/collectSampleMetrics.log"
    errorCode=$?

    # Verify collectSampleMetrics.sh error handling behavior
    assertEquals "collectSampleMetrics.sh returned incorrect error code when the SAM file was corrupt." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running collectSampleMetrics.sh."
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Error detected while running collectSampleMetrics.sh."
    assertFileContains "$tempDir/error.log" "samtools view"
    assertFileContains "$logDir/collectSampleMetrics.log" "sam file"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "collectSampleMetrics.sh finished"
    assertFileNotContains "$logDir/collectSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the collectSampleMetrics.sh script traps errors
testCollectSampleMetricsSamtoolsViewTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCollectSampleMetricsSamtoolsViewTrap 100
}

# Verify the collectSampleMetrics.sh script traps errors
testCollectSampleMetricsSamtoolsViewTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCollectSampleMetricsSamtoolsViewTrap 98
}

# Verify the collectSampleMetrics.sh script traps errors
testCollectSampleMetricsSamtoolsViewTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCollectSampleMetricsSamtoolsViewTrap 100
}


# Verify the combineSampleMetrics.sh script detects missing input file
tryCombineSampleMetricsMissingSampleDirRaiseGlobalError()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Run create_snp_list.py with missing sampleDirectories.txt
    combineSampleMetrics.sh "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "combineSampleMetrics.sh returned incorrect error code when sample directories file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "combineSampleMetrics.sh failed."
    assertFileNotContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh failed."
    assertFileContains "$tempDir/error.log" "Sample directories file $tempDir/sampleDirectories.txt does not exist"
    assertFileContains "$logDir/combineSampleMetrics.log" "Sample directories file $tempDir/sampleDirectories.txt does not exist"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "create_snp_list.py finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the combineSampleMetrics.sh script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the combineSampleMetrics.sh script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}

# Verify the combineSampleMetrics.sh script detects missing input file
testCombineSampleMetricsMissingSampleDirRaiseGlobalErrorStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsMissingSampleDirRaiseGlobalError 100
}


# Verify the combineSampleMetrics.sh script detects a missing sample metrics file
tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Setup directories and env variables used to trigger error handling.
    # This simulates what run_snp_pipeline does before running other scripts
    export logDir="$tempDir/logs"
    mkdir -p "$logDir"
    export errorOutputFile="$tempDir/error.log"

    # Try to combine metrics
    printf "%s\n" $tempDir/samples/* > "$tempDir/sampleDirectories.txt"
    echo "Dummy snpma.fasta content" > "$tempDir/snpma.fasta"
    combineSampleMetrics.sh -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify combineSampleMetrics.sh error handling behavior
    assertEquals "combineSampleMetrics.sh returned incorrect error code when the sample metrics file was missing." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "combineSampleMetrics.sh warning"
    assertFileNotContains "$tempDir/error.log" "combineSampleMetrics.sh failed"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh failed"
    assertFileContains "$tempDir/error.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$logDir/combineSampleMetrics.log" "Sample metrics file $tempDir/samples/sample1/metrics does not exist"
    assertFileContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
}

# Verify the combineSampleMetrics.sh script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the combineSampleMetrics.sh script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}

# Verify the combineSampleMetrics.sh script detects a missing sample metrics file
testCombineSampleMetricsMissingSampleMetricsRaiseSampleWarningStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsMissingSampleMetricsRaiseSampleWarning 0
}


# Verify the combineSampleMetrics.sh script traps attempts to write to unwritable file
tryCombineSampleMetricsPermissionTrap()
{
    expectErrorCode=$1
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    combineSampleMetrics.sh -o "$tempDir/metrics.tsv" "$tempDir/sampleDirectories.txt" &> "$logDir/combineSampleMetrics.log"
    errorCode=$?

    # Verify combineSampleMetrics.sh error handling behavior
    assertEquals "combineSampleMetrics.sh returned incorrect error code when the output file was unwritable." $expectErrorCode $errorCode
    verifyNonEmptyReadableFile "$tempDir/error.log"
    assertFileContains "$tempDir/error.log" "Error detected while running combineSampleMetrics.sh"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Error detected while running combineSampleMetrics.sh"
    assertFileContains "$tempDir/error.log" "printf"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"
    assertFileNotContains "$logDir/combineSampleMetrics.log" "Use the -f option to force a rebuild"
    rm -f "$tempDir/metrics.tsv"
}

# Verify the combineSampleMetrics.sh script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapStop()
{
    export SnpPipeline_StopOnSampleError=true
    tryCombineSampleMetricsPermissionTrap 100
}

# Verify the combineSampleMetrics.sh script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapNoStop()
{
    export SnpPipeline_StopOnSampleError=false
    tryCombineSampleMetricsPermissionTrap 100
}

# Verify the combineSampleMetrics.sh script traps attempts to write to unwritable file
testCombineSampleMetricsPermissionTrapStopUnset()
{
    unset SnpPipeline_StopOnSampleError
    tryCombineSampleMetricsPermissionTrap 100
}


# Verify the run_snp_pipeline.sh script detects missing reference
tryRunSnpPipelineMissingReferenceRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}

# Verify the collectSampleMetrics.sh script detects missing reference
testRunSnpPipelineMissingReferenceRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingReferenceRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing configuration file
tryRunSnpPipelineMissingConfigurationFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing configuration file
testRunSnpPipelineMissingConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingConfigurationFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir
    copy_snppipeline_data.py configurationFile $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects invalid aligner configuration
testRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineInvalidAlignerConfigurationFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing file of sample directories
tryRunSnpPipelineMissingSampleDirFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing file of sample directories
testRunSnpPipelineMissingSampleDirFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSampleDirFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script validates the file of sample directories
tryRunSnpPipelineValidateSampleDirFileRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineValidateSampleDirFileRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir
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
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
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

    assertFileContains "$tempDir/snpma.fasta" "sample4"
    assertFileContains "$tempDir/snpma.fasta" "sample5"
    assertFileContains "$tempDir/snpma.vcf" "sample4"
    assertFileContains "$tempDir/snpma.vcf" "sample5"

    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "See the log file $tempDir/error.log for a summary of errors."
}

# Verify the run_snp_pipeline.sh script validates the file of sample directories
testRunSnpPipelineValidateSampleDirFileRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineValidateSampleDirFileRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects missing directory of samples
tryRunSnpPipelineMissingSamplesDirRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects missing directory of samples
testRunSnpPipelineMissingSamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineMissingSamplesDirRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects empty directory of samples
tryRunSnpPipelineEmptySamplesDirRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects empty directory of samples
testRunSnpPipelineEmptySamplesDirRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineEmptySamplesDirRaiseFatalError 1
}


# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError()
{
    expectErrorCode=$1

    # Extract test data to temp dir
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}

# Verify the run_snp_pipeline.sh script detects directory of samples with no fastq files in any subdirectories
testRunSnpPipelineSamplesDirNoFastqRaiseFatalErrorStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineSamplesDirNoFastqRaiseFatalError 1
}


# Verify run_snp_pipeline.sh trap handling
tryRunSnpPipelineTrapPrepReferenceTrap()
{
    expectErrorCode=$1

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Deliberately corrupt the fasta file
    sed -i 's/>/@@@/g' "$tempDir/reference/lambda_virus.fasta"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when the input fasta was corrupt." $expectErrorCode $errorCode
    assertFileContains "$tempDir/error.log" "Error detected while running prepReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2-build"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running prepReference.sh."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2-build"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepReference.sh finished"
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
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.snp.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapNoStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapPrepReferenceTrapStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapPrepReferenceTrap 1
}


# Verify run_snp_pipeline.sh trap handling
tryRunSnpPipelineTrapAlignSampleToReferenceTrap()
{
    expectErrorCode=$1

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "prepSamples.sh"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "create_snp_list.py"

    assertFileContains "$tempDir/error.log" "Shutting down the SNP Pipeline"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Shutting down the SNP Pipeline"
    assertFileNotContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples"

    verifyNonExistingFile "$tempDir/samples/sample1/reads.unsorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.sorted.bam"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.snp.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapStop()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=true" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapAlignSampleToReferenceTrap 1
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapNoStopAllFail()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    expectErrorCode=1

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "bowtie2"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "prepSamples.sh"

    assertFileContains "$tempDir/error.log" "prepSamples.sh failed"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample2/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample3/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample4/reads.sam"

    assertFileContains "$tempDir/error.log" "create_snp_list.py failed"
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
    verifyNonExistingFile "$tempDir/samples/sample1/reads.all.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/var.flt.vcf"
    verifyNonExistingFile "$tempDir/snplist.txt"
    verifyNonExistingFile "$tempDir/samples/sample1/reads.snp.pileup"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.fasta"
    verifyNonExistingFile "$tempDir/samples/sample1/consensus.vcf"
    verifyNonExistingFile "$tempDir/snpma.fasta"
    verifyNonExistingFile "$tempDir/snpma.vcf"
    verifyNonExistingFile "$tempDir/referenceSNP.fasta"
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapNoStopSomeFail()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"
    expectErrorCode=0

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Deliberately corrupt some of the fastq files
    echo "Garbage" > "$tempDir/samples/sample1/sample1_1.fastq"
    echo "Garbage" > "$tempDir/samples/sample4/sample4_1.fastq"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" 2> "$tempDir/run_snp_pipeline.stderr.log" > "$tempDir/run_snp_pipeline.stdout.log"
    errorCode=$?

    # Verify error handling behavior
    assertEquals "run_snp_pipeline.sh did not return an error code when some the input fastq files were corrupt." $expectErrorCode $errorCode
    assertFileContains "$tempDir/error.log" "Error detected while running alignSampleToReference.sh."
    assertFileContains "$tempDir/error.log" "bowtie2"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample1/sample1_1.fastq $tempDir/samples/sample1/sample1_2.fastq"
    assertFileNotContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample2/sample2_1.fastq $tempDir/samples/sample2/sample2_2.fastq"
    assertFileNotContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample3/sample3_1.fastq $tempDir/samples/sample3/sample3_2.fastq"
    assertFileContains "$tempDir/error.log" "alignSampleToReference.sh $tempDir/reference/lambda_virus.fasta $tempDir/samples/sample4/sample4_1.fastq $tempDir/samples/sample4/sample4_2.fastq"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "alignSampleToReference.sh finished"
    assertFileNotContains "$tempDir/run_snp_pipeline.stdout.log" "Use the -f option to force a rebuild"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "prepSamples.sh"

    assertFileContains "$tempDir/error.log" "prepSamples.sh failed"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample1/reads.sam"
    assertFileNotContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample2/reads.sam"
    assertFileNotContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample3/reads.sam"
    assertFileContains "$tempDir/error.log" "Sample SAM file $tempDir/samples/sample4/reads.sam"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "prepSamples.sh finished"

    assertFileNotContains "$tempDir/error.log" "create_snp_list.py failed"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample1/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "VCF file $tempDir/samples/sample4/var.flt.vcf does not exist"
    assertFileContains "$tempDir/error.log" "Error: 2 VCF files were missing or empty"
    assertFileContains "$tempDir/run_snp_pipeline.stdout.log" "create_snp_list.py finished"

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

    assertFileContains "$tempDir/snpma.fasta" "sample2"
    assertFileContains "$tempDir/snpma.fasta" "sample3"
    assertFileContains "$tempDir/snpma.vcf" "sample2"
    assertFileContains "$tempDir/snpma.vcf" "sample3"

    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "There were errors processing some samples."
    assertFileContains "$tempDir/run_snp_pipeline.stderr.log" "See the log file $tempDir/error.log for a summary of errors."
}

# Verify run_snp_pipeline.sh trap handling
testRunSnpPipelineTrapAlignSampleToReferenceTrapStopUnset()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")
    copy_snppipeline_data.py configurationFile $tempDir
    echo "unset SnpPipeline_StopOnSampleError" >> "$tempDir/snppipeline.conf"
    tryRunSnpPipelineTrapAlignSampleToReferenceTrap 1
}


# Verify run_snp_pipeline generates correct results for the lambda data set
testRunSnpPipelineLambda()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir/originalInputs

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
    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snplist.txt"                   "$tempDir/expectedResults/snplist.txt"
    assertIdenticalFiles "$tempDir/snpma.fasta"                   "$tempDir/expectedResults/snpma.fasta"
    assertIdenticalFiles "$tempDir/snpma.vcf"                     "$tempDir/expectedResults/snpma.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source --ignore-matching-lines=##bcftools
    assertIdenticalFiles "$tempDir/referenceSNP.fasta"            "$tempDir/expectedResults/referenceSNP.fasta"
    assertIdenticalFiles "$tempDir/metrics.tsv"                   "$tempDir/expectedResults/metrics.tsv"
    assertIdenticalFiles "$tempDir/samples/sample1/consensus.vcf" "$tempDir/expectedResults/samples/sample1/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample2/consensus.vcf" "$tempDir/expectedResults/samples/sample2/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample3/consensus.vcf" "$tempDir/expectedResults/samples/sample3/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source
    assertIdenticalFiles "$tempDir/samples/sample4/consensus.vcf" "$tempDir/expectedResults/samples/sample4/consensus.vcf" --ignore-matching-lines=##fileDate --ignore-matching-lines=##source

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
}


# Verify run_snp_pipeline generates correct results for the lambda data set
testRunSnpPipelineLambdaUnpaired()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir/originalInputs

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
    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-1" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-2" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-3" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-4" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/prepSamples.log-1" "prepSamples.sh finished"
    assertFileContains "$logDir/prepSamples.log-2" "prepSamples.sh finished"
    assertFileContains "$logDir/prepSamples.log-3" "prepSamples.sh finished"
    assertFileContains "$logDir/prepSamples.log-4" "prepSamples.sh finished"
    assertFileContains "$logDir/snpList.log" "create_snp_list.py finished"
    assertFileContains "$logDir/callConsensus.log-1" "call_consensus.py finished"
    assertFileContains "$logDir/callConsensus.log-2" "call_consensus.py finished"
    assertFileContains "$logDir/callConsensus.log-3" "call_consensus.py finished"
    assertFileContains "$logDir/callConsensus.log-4" "call_consensus.py finished"
    assertFileContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileContains "$logDir/snpMatrix.log" "create_snp_matrix.py finished"
    assertFileContains "$logDir/snpReference.log" "create_snp_reference_seq.py finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/collectSampleMetrics.log-2" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/collectSampleMetrics.log-4" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/collectSampleMetrics.log-3" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"
}


# Verify run_snp_pipeline runs to completion with a single sample
testRunSnpPipelineLambdaSingleSample()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir/originalInputs

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
    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-1" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/prepSamples.log-1" "prepSamples.sh finished"
    assertFileContains "$logDir/snpList.log" "create_snp_list.py finished"
    assertFileContains "$logDir/callConsensus.log-1" "call_consensus.py finished"
    assertFileContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileContains "$logDir/snpMatrix.log" "create_snp_matrix.py finished"
    assertFileContains "$logDir/snpReference.log" "create_snp_reference_seq.py finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"

    # Verify correct results
    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snpma.vcf"                     "$tempDir/samples/sample1/consensus.vcf"  # Just copy the sample VCF to the snpma.vcf
}


# Verify run_snp_pipeline runs to completion with no snps and works properly when re-run
testRunSnpPipelineZeroSnps()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    touch -d  '-4 day' $tempDir/samples/*/reads.all.pileup
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

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/prepReference.log" "prepReference.sh finished"
    assertFileContains "$logDir/alignSamples.log-1" "alignSampleToReference.sh finished"
    assertFileContains "$logDir/prepSamples.log-1" "prepSamples.sh finished"
    assertFileContains "$logDir/snpList.log" "create_snp_list.py finished"
    assertFileContains "$logDir/callConsensus.log-1" "call_consensus.py finished"
    assertFileContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileContains "$logDir/snpMatrix.log" "create_snp_matrix.py finished"
    assertFileContains "$logDir/snpReference.log" "create_snp_reference_seq.py finished"
    assertFileContains "$logDir/collectSampleMetrics.log-1" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"
}


# Verify run_snp_pipeline rebuilds the snplist when at least one var.flt.vcf is missing and 
# at least one var.flt.vcf is newer
testRunSnpPipelineRerunMissingVCF()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

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
    touch -d  '-4 day' $tempDir/samples/*/reads.all.pileup
    touch $tempDir/samples/*/var.flt.vcf
    rm -rf "$tempDir/samples/sample1"
    sleep 1

    copy_snppipeline_data.py configurationFile $tempDir
    echo "SnpPipeline_StopOnSampleError=false" >> "$tempDir/snppipeline.conf"

    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -S "$tempDir/sampleDirectories.txt" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify output results
    verifyNonEmptyReadableFile "$tempDir/snplist.txt"
    assertNewerFile "$tempDir/snplist.txt" "$tempDir/samples/sample2/var.flt.vcf"
    assertFileNotContains "$tempDir/snplist.txt" "sample1"

    # Verify output results exist, and no snps were found
    verifyNonEmptyReadableFile "$tempDir/snpma.fasta"
    verifyNonEmptyReadableFile "$tempDir/snpma.vcf"
    verifyNonEmptyReadableFile "$tempDir/referenceSNP.fasta"

    # Verify log files
    logDir=$(echo $(ls -d $tempDir/logs*))
    verifyNonEmptyReadableFile "$logDir/snppipeline.conf"
    assertFileContains "$logDir/prepSamples.log-2" "prepSamples.sh finished"
    assertFileContains "$logDir/snpList.log" "create_snp_list.py finished"
    assertFileContains "$logDir/callConsensus.log-2" "call_consensus.py finished"
    assertFileContains "$logDir/mergeVcf.log" "mergeVcf.sh finished"
    assertFileContains "$logDir/snpMatrix.log" "create_snp_matrix.py finished"
    assertFileContains "$logDir/snpReference.log" "create_snp_reference_seq.py finished"
    assertFileContains "$logDir/collectSampleMetrics.log-2" "collectSampleMetrics.sh finished"
    assertFileContains "$logDir/combineSampleMetrics.log" "combineSampleMetrics.sh finished"
}


# Verify processing steps are skipped when output files are already fresh.
testAlreadyFreshOutputs()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> /dev/null

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Force timestamps to change so outputs are newer than inputs.
    # The files are small, quickly processed, and timestamps might not differ when we expect they will differ.
    touch -d '-10 day' $tempDir/reference/*.fasta
    touch -d  '-9 day' $tempDir/reference/*.bt2
    touch -d  '-8 day' $tempDir/samples/*/*.fastq
    touch -d  '-7 day' $tempDir/samples/*/reads.sam
    touch -d  '-6 day' $tempDir/samples/*/reads.unsorted.bam
    touch -d  '-5 day' $tempDir/samples/*/reads.sorted.bam
    touch -d  '-4 day' $tempDir/samples/*/reads.all.pileup
    touch -d  '-3 day' $tempDir/samples/*/var.flt.vcf
    touch -d  '-2 day' $tempDir/snplist.txt
    touch -d  '-1 day' $tempDir/samples/*/consensus.vcf

    # Test special collectSampleMetrics result persistence
    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=20000"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=94.54"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=22.89"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=286.84"
    echo numberReads=AAAAA > "$tempDir/samples/sample1/metrics"
    echo percentReadsMapped=BBBBB >> "$tempDir/samples/sample1/metrics"
    echo avePileupDepth=CCCCC >> "$tempDir/samples/sample1/metrics"
    echo aveInsertSize=DDDDD >> "$tempDir/samples/sample1/metrics"

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
    assertFileContains "$logDir/prepSamples.log-1" "Vcf file is already freshly created for sample1.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-2" "Vcf file is already freshly created for sample2.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-3" "Vcf file is already freshly created for sample4.  Use the -f option to force a rebuild."
    assertFileContains "$logDir/prepSamples.log-4" "Vcf file is already freshly created for sample3.  Use the -f option to force a rebuild."

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

    # Special collectSampleMetrics re-use last metrics
    assertFileNotContains "$tempDir/samples/sample1/metrics" "numberReads=20000"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=94.54"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "avePileupDepth=22.89"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "aveInsertSize=286.84"
    assertFileContains "$tempDir/samples/sample1/metrics" "numberReads=AAAAA"
    assertFileContains "$tempDir/samples/sample1/metrics" "percentReadsMapped=BBBBB"
    assertFileContains "$tempDir/samples/sample1/metrics" "avePileupDepth=CCCCC"
    assertFileContains "$tempDir/samples/sample1/metrics" "aveInsertSize=DDDDD"
    assertFileContains "$tempDir/metrics.tsv" "sample1.*AAAAA.*BBBBB.*DDDDD.*CCCCC"
}


# Verify underscores in metrics.tsv column headers can be controlled with configuration file
testRunSnpPipelineMetricsColumnHeadingsUnderscores()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir/originalInputs

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -m copy -o "$tempDir" -s "$tempDir/originalInputs/samples" "$tempDir/originalInputs/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify output metrics
    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/metrics.tsv" "$tempDir/expectedResults/metrics.tsv"
    head -n 1 "$tempDir/metrics.tsv" | grep "_" > /dev/null
    assertTrue "No underscores were found in the metrics column headings"  $?

    # Delete the metrics file and re-run with the option to use spaces
    rm "$tempDir/metrics.tsv"
    copy_snppipeline_data.py configurationFile $tempDir
    echo 'CombineSampleMetrics_ExtraParams="-s"' >> "$tempDir/snppipeline.conf"
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify output metrics have no underscores
    head -n 1 "$tempDir/metrics.tsv" | grep "_" > /dev/null
    assertFalse "Underscores should not be found in the metrics column headings when using -s combineSampleMetrics option"  $?
}


# Verify samples with excessive snps are excluded from the snplist, snp matrix, and snpma.vcf.
testRunSnpPipelineExcessiveSnps()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py lambdaVirusInputs $tempDir

    # Create a config file with a low enough maxsnps setting to block a sample
    copy_snppipeline_data.py configurationFile $tempDir
    sed -i s:SnpPipeline_MaxSnps=-1:SnpPipeline_MaxSnps=45: "$tempDir/snppipeline.conf"

    # Run the pipeline, specifing the locations of samples and the reference
    run_snp_pipeline.sh -c "$tempDir/snppipeline.conf" -o "$tempDir" -s "$tempDir/samples" "$tempDir/reference/lambda_virus.fasta" &> "$tempDir/run_snp_pipeline.log"

    # Verify no errors
    verifyNonExistingFile "$tempDir/error.log"

    # Verify no weird freshness skips
    assertFileNotContains "$tempDir/run_snp_pipeline.log" "Use the -f option to force a rebuild"

    # Verify output
    assertFileContains "$tempDir/samples/sample1/metrics" "excludedSample=Excluded$"
    assertFileContains "$tempDir/samples/sample2/metrics" "excludedSample=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "snps=$"
    assertFileContains "$tempDir/samples/sample2/metrics" "snps=44$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "missingPos=$"
    assertFileContains "$tempDir/samples/sample1/metrics" "errorList=.*Excluded: exceeded 45 maxsnps."
    assertFileContains "$tempDir/samples/sample2/metrics" "errorList=\"No compressed fastq.gz or fq.gz files were found.\"$"
    assertFileContains "$tempDir/metrics.tsv"             "sample1.*Excluded.*Excluded: exceeded 45 maxsnps."
    assertFileNotContains "$tempDir/samples/sample1/metrics" "Consensus.*not found"
    assertFileNotContains "$tempDir/samples/sample1/metrics" "Consensus.*not found"
    assertFileNotContains "$tempDir/metrics.tsv"             "sample1.*Consensus.*not found"
    assertFileNotContains "$tempDir/snplist.txt" "sample1"
    assertFileNotContains "$tempDir/snpma.fasta" "sample1"
    assertFileNotContains "$tempDir/snpma.vcf"   "sample1"

    copy_snppipeline_data.py lambdaVirusExpectedResults $tempDir/expectedResults
    grep -v sample1 "$tempDir/expectedResults/metrics.tsv" > "$tempDir/expectedResults/metrics.withoutSample1.tsv"
    grep -v sample1 "$tempDir/metrics.tsv" > "$tempDir/metrics.withoutSample1.tsv"
    assertIdenticalFiles "$tempDir/metrics.withoutSample1.tsv" "$tempDir/expectedResults/metrics.withoutSample1.tsv"
}


# Verify run_snp_pipeline generates correct results for the Salmonella Agona data set
testRunSnpPipelineAgona()
{
    tempDir=$(mktemp -d -p "$SHUNIT_TMPDIR")

    # Copy the supplied test data to a work area:
    copy_snppipeline_data.py agonaInputs "$tempDir"

    # Create sample directories
    mkdir -p "$tempDir/samples/ERR178926"  "$tempDir/samples/ERR178927"  "$tempDir/samples/ERR178928"  "$tempDir/samples/ERR178929"  "$tempDir/samples/ERR178930"

    # Download sample data from SRA at NCBI.
    fastq-dump --split-files --outdir "$tempDir/samples/ERR178926" ERR178926 &> /dev/null
    fastq-dump --split-files --outdir "$tempDir/samples/ERR178927" ERR178927 &> /dev/null
    fastq-dump --split-files --outdir "$tempDir/samples/ERR178928" ERR178928 &> /dev/null
    fastq-dump --split-files --outdir "$tempDir/samples/ERR178929" ERR178929 &> /dev/null
    fastq-dump --split-files --outdir "$tempDir/samples/ERR178930" ERR178930 &> /dev/null

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
    copy_snppipeline_data.py agonaExpectedResults $tempDir/expectedResults
    assertIdenticalFiles "$tempDir/snplist.txt"        "$tempDir/expectedResults/snplist.txt"
    assertIdenticalFiles "$tempDir/snpma.fasta"        "$tempDir/expectedResults/snpma.fasta"
    assertIdenticalFiles "$tempDir/referenceSNP.fasta" "$tempDir/expectedResults/referenceSNP.fasta"
    assertIdenticalFiles "$tempDir/metrics.tsv"        "$tempDir/expectedResults/metrics.tsv"
}


# load shunit2 and execute all the tests in this script
. shunit2
