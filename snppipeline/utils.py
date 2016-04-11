from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pprint
import sys
import time
import traceback
import vcf


#==============================================================================
#Prep work
#==============================================================================

verbose_print  = lambda *a, **k: None
verbose_pprint = lambda *a, **k: None


def set_logging_verbosity(options_dict):
    """Enable or disable logging.

    Args:
        verbose : Verbosity value, any value greater than 0 enables logging
    """
    global verbose_print
    global verbose_pprint
    verbose_print  = print         if options_dict['verbose'] > 0 else lambda *a, **k: None
    verbose_pprint = pprint.pprint if options_dict['verbose'] > 0 else lambda *a, **k: None


#==============================================================================
#Define functions
#==============================================================================

def timestamp():
    """Return a timestamp string."""
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())


def program_name():
    """Return the basename of the python script being executed."""
    return os.path.basename(sys.argv[0])


def command_line_short():
    """Return the command line string without the full path to the program."""
    return "%s %s" % (program_name(), " ".join(sys.argv[1:]))


def command_line_long():
    """Return the command line string with the full path to the program."""
    return " ".join(sys.argv)


def global_error(message):
    """
    Log a fatal error to the error summary file and exit with error code 100
    to cause Sun Grid Engine to also detect the error.

    Args:
        message : str
            Error message
    """
    # Log the event to the error log
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            print("%s failed." % program_name(), file=err_log)
            if message:
                print(message, file=err_log)
            print("=" * 80, file=err_log)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    if message:
        print(message, file=sys.stderr)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    sys.exit(100)


def sample_error(message, continue_possible=False):
    """
    Log a fatal error to the error summary file and exit with error code 100
    to cause Sun Grid Engine to also detect the error.

    Args:
        message : str
            Error message
        continue_possible : boolean
            Indicates if it is possible to continue execution.  Setting this
            flag true may allow the code to continue withou exiting if
            configured to do so.
    """
    stop_on_error_env = os.environ.get("SnpPipeline_StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"

    # Log the event to the error log
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            if stop_on_error or not continue_possible:
                print("%s failed." % program_name(), file=err_log)
            else:
                print("%s" % program_name(), file=err_log)
            print(message, file=err_log)
            print("=" * 80, file=err_log)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    print(message, file=sys.stderr)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    if stop_on_error:
        # run_snp_pipeline.sh will know this error has already been reported
        sys.exit(100)
    else:
        # run_snp_pipeline.sh will know this error has already been reported,
        # but it should not stop execution
        if not continue_possible:
            sys.exit(98)


def handle_global_exception(exc_type, exc_value, exc_traceback):
    """
    This function replaces the default python unhandled exception handler.
    It Logs the error and returns error code 100 to cause Sun Grid Engine to
    also detect the error.
    """
    # Report the exception in the error log if configuired
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            trace_entries = traceback.extract_tb(exc_traceback)
            file_name, line_number, function_name, code_text = trace_entries[-1]
            exc_type_name = exc_type.__name__

            print("Error detected while running %s." % program_name(), file=err_log)
            print("", file=err_log)
            print("The command line was:", file=err_log)
            print("    %s" % command_line_short(), file=err_log)
            print("", file=err_log)
            print("%s exception in function %s at line %d in file %s" % (exc_type_name, function_name, line_number, file_name), file=err_log)
            print("    %s" % code_text, file=err_log)
            print("=" * 80, file=err_log)

    # Report the exception in the usual way to stderr
    sys.stdout.flush() # make sure stdout is flushed before printing the trace
    traceback.print_exception(exc_type, exc_value, exc_traceback)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    sys.exit(100)


def handle_sample_exception(exc_type, exc_value, exc_traceback):
    """
    This function replaces the default python unhandled exception handler.
    It Logs the error and returns error code 100 to cause Sun Grid Engine to
    also detect the error.
    """
    # Report the exception in the error log if configuired
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            trace_entries = traceback.extract_tb(exc_traceback)
            file_name, line_number, function_name, code_text = trace_entries[-1]
            exc_type_name = exc_type.__name__

            print("Error detected while running %s." % program_name(), file=err_log)
            print("", file=err_log)
            print("The command line was:", file=err_log)
            print("    %s" % command_line_short(), file=err_log)
            print("", file=err_log)
            print("%s exception in function %s at line %d in file %s" % (exc_type_name, function_name, line_number, file_name), file=err_log)
            print("    %s" % code_text, file=err_log)
            print("=" * 80, file=err_log)

    # Report the exception in the usual way to stderr
    sys.stdout.flush() # make sure stdout is flushed before printing the trace
    traceback.print_exception(exc_type, exc_value, exc_traceback)

    # Exit 100 does two things:
    # 1. Sun Grid Engine will stop execution of dependent jobs
    # 2. run_snp_pipeline.sh will know this error has already been reported
    stop_on_error_env = os.environ.get("SnpPipeline_StopOnSampleError")
    stop_on_error = stop_on_error_env is None or stop_on_error_env == "true"
    if stop_on_error:
        # run_snp_pipeline.sh will know this error has already been reported
        sys.exit(100)
    else:
        # run_snp_pipeline.sh will know this error has already been reported,
        # but it should not stop execution
        sys.exit(98)


def report_error(message):
    """Send an error message to the error log and to stderr.
    """
    error_output_file = os.environ.get("errorOutputFile")
    if error_output_file:
        with open(error_output_file, "a") as err_log:
            print(message, file=err_log)

    # send the detail error message to stderr -- this will put the error
    # message in the process-specific log file.
    sys.stdout.flush() # make sure stdout is flushed before printing the error
    if message:
        print(message, file=sys.stderr)


def verify_existing_input_files(error_prefix, file_list):
    """Verify each file in a list of files exists.  It does
    not matter whether the file is empty.
    Missing files are reported in the verbose log.

    Args:
        error_prefix : first part of error message to be logged
        file_list : list of relative or absolute paths to files

    Returns:
        int number of missing files
    """
    bad_count = 0
    for file_path in file_list:

        if not os.path.isfile(file_path):
            bad_count += 1
            err_message = "%s %s does not exist." % (error_prefix, file_path)
            report_error(err_message)
            continue

    return bad_count


def verify_non_empty_input_files(error_prefix, file_list):
    """Verify each file in a list of files exists and is non-empty.
    Missing or empty files are reported in the verbose log.

    Args:
        error_prefix : first part of error message to be logged
        file_list : list of relative or absolute paths to files

    Returns:
        int number of missing or empty files
    """
    bad_count = 0
    for file_path in file_list:

        if not os.path.isfile(file_path):
            bad_count += 1
            err_message = "%s %s does not exist." % (error_prefix, file_path)
            report_error(err_message)
            continue
        if os.path.getsize(file_path) == 0:
            bad_count += 1
            err_message = "%s %s is empty." % (error_prefix, file_path)
            report_error(err_message)
            continue

    return bad_count


def target_needs_rebuild(source_files, target_file):
    """Determine if a target file needs a fresh rebuild, i.e. the target does
    not exist or its modification time is older than any of its source files.

    Args:
        source_files : relative or absolute path to a list of files
        target_file : relative or absolute path to target file
    """
    if not os.path.isfile(target_file):
        return True

    if os.path.getsize(target_file) == 0:
        return True

    target_timestamp = os.stat(target_file).st_mtime

    for source_file in source_files:
        # A non-existing source file should neither force a rebuild, nor prevent a rebuild.
        # You should error-check the existence of the source files before calling this function.
        #
        # An empty source file should force a rebuild if it is newer than the target, just like
        # a regular non-empty source file.
        if not os.path.isfile(source_file):
            continue

        source_timestamp = os.stat(source_file).st_mtime
        if source_timestamp > target_timestamp:
            return True

    return False


def create_snp_pileup(all_pileup_file_path, snp_pileup_file_path, snp_set):
    """Create a subset pileup with SNP locations only.

    Given a whole-genome pileup file, create a new pileup file with a subset
    of the records at the SNP locations only.

    Args:
        all_pileup_file_path: path to a whole-genome pileup file
        all_pileup_file_path: path to a snp pileup file to be created
        snp_set: set of (CHROM, POS) tuples identifying the locations with SNPs
    """
    with open(all_pileup_file_path, "r") as all_pileup_file_object:
        with open(snp_pileup_file_path, "w") as snp_pileup_file_object:
            for pileup_line in all_pileup_file_object:
                current_line_data = pileup_line.rstrip().split()
                seq_id, pos = current_line_data[:2]
                key = (seq_id, int(pos))
                if key in snp_set:
                    snp_pileup_file_object.write(pileup_line)


def write_list_of_snps(file_path, snp_dict):
    """Write out list of snps for all samples to a single file.

    Args:
        file_path : path to snplist file to be written
        snp_dict  : dictionary with key = tuple(CHROM, POS) -> value = list[sampleName1, sampleName2, ..., sampleNameN]

    Returns:
        Nothing
    """

    with open(file_path, "w") as snp_list_file_object:
        for key in sorted(snp_dict.keys()):
            sample_list = snp_dict[key]
            snp_list_file_object.write("%s\t%d\t%d\t%s\n" % (key[0], key[1], len(sample_list), "\t".join(sample_list)))


def read_snp_position_list(snp_list_file_path):
    """Read list of snp positions across all samples from the snplist.txt.

    Args:
        snp_list_file_path : path to snplist file to be written

    Returns:
        snp_list  : sorted list of tuple(str(CHROM), int(POS))
    """

    snp_list = list()
    with open(snp_list_file_path, "r") as snp_list_file_object:
        for line in snp_list_file_object:
            chrom, pos = line.split()[0:2]
            snp_list.append((chrom, int(pos)))
    return snp_list


def write_reference_snp_file(reference_file_path, snp_list_file_path,
                             snp_reference_file_path):
    """Write out the snp fasta file for the reference.fasta using the snp
    position file ( snplist.txt).
    """
    #TODO finish documentation
    #TODO actual code is more general than stated. Fix this.

    with open(snp_list_file_path, "r") as snp_list_file:
        position_list = [line.split()[0:2] for line in snp_list_file]
    match_dict = SeqIO.to_dict(SeqIO.parse(reference_file_path, "fasta"))

    with open(snp_reference_file_path, "w") as snp_reference_file_object:
        for ordered_id in sorted(match_dict.keys()):
            ref_str = ""
            for chrom_id, pos in position_list:
                if chrom_id == ordered_id:
                    ref_str += match_dict[ordered_id][int(pos) - 1].upper()
            record = SeqRecord(Seq(ref_str), id=ordered_id, description="")
            SeqIO.write([record], snp_reference_file_object, "fasta")


def convert_vcf_file_to_snp_set(vcf_file_path):
    """convert vcf files to a set of SNPs.

    Args:
        vcf_file_path : relative or absolute path to the sample VCF file

    Returns:
        snp_set  : set of (CHROM, POS) tuples

    """

    snp_set = set()

    with open(vcf_file_path, 'r') as vcf_file_object:
        vcf_reader = vcf.Reader(vcf_file_object)
        for vcf_data_line in vcf_reader:
            key = (vcf_data_line.CHROM, vcf_data_line.POS)
            snp_set.add(key)

    return snp_set


def calculate_sequence_distance(seq1, seq2, case_insensitive=True):
    """Calulate the number of nucleotide differences between two sequences.

    The sequences must be the same length.

    Args:
        seq1 : DNA string 1
        seq2 : DNA string 2
        case_insensitive : optional flag for case insensitive compare, defaults to True

    Returns:
        int number of differences
    """
    if case_insensitive:
        allowed_bases = frozenset(['A', 'C', 'G', 'T'])
        seq1 = seq1.upper()
        seq2 = seq2.upper()
    else:
        allowed_bases = frozenset(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'])

    mismatches = 0
    for pos in range(len(seq1)):
        base1 = seq1[pos]
        base2 = seq2[pos]
        if base1 not in allowed_bases:
            continue
        if base2 not in allowed_bases:
            continue
        if base1 != base2:
            mismatches += 1
    return mismatches
