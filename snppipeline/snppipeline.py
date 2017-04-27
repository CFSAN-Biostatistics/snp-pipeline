from __future__ import print_function
from __future__ import absolute_import
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import os
from snppipeline import pileup
from snppipeline import utils
from snppipeline import vcf_writer


def create_snp_list(args):
    """Create SNP list file

    Description:
    Create the SNP list -- the list of positions where variants were found
    and the corresponding list of samples having a variant at each position.
    This function expects, or creates '(*)', the following files arranged
    in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/var.flt.vcf
                ...
            snplist.txt (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to
           the sample directories.
        2. The var.flt.vcf variant input files are used to construct the
           SNP position list.
        3. The snplist.txt output file contains the union of the SNP positions
           and sample names extracted from all the var.flt.vcf files.

    The sampleDirectories.txt and var.flt.vcf files are created outside of
    this function. The package documentation provides an example of creating
    these files based on the lambda_virus sequence that is used as one test
    for this package.

    Args:
        sampleDirsFile: File path (not just file name) of file containing paths
            to directories containing var.flt.vcf file for each sequence.
        vcfFileName: File name of the VCF files which must exist in each of the
            sample directories
        snpListFile: File path (not just file name) of text format list
            of SNP positions

    Raises:

    Examples:
    args = argparse.Namespace
    args.sampleDirsFile = 'sampleDirectories.txt'
    args.vcfFileName = 'var.flt.vcf'
    args.snpListFile = 'snplist.txt'
    create_snp_list(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Prep work
    #==========================================================================
    sample_directories_list_path = args.sampleDirsFile
    bad_file_count = utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_path])
    if bad_file_count > 0:
        utils.global_error(None)

    with open(sample_directories_list_path, "r") as sample_directories_list_file:
        unsorted_list_of_sample_directories = [line.rstrip() for line in sample_directories_list_file]
    unsorted_list_of_sample_directories = [d for d in unsorted_list_of_sample_directories if d]
    sorted_list_of_sample_directories = sorted(unsorted_list_of_sample_directories)

    #==========================================================================
    # Validate inputs
    #==========================================================================
    snp_list_file_path = args.snpListFile
    vcf_file_name = args.vcfFileName
    list_of_vcf_files = [os.path.join(dir, vcf_file_name) for dir in sorted_list_of_sample_directories]

    bad_file_count = utils.verify_non_empty_input_files("VCF file", list_of_vcf_files)
    if bad_file_count == len(list_of_vcf_files):
        utils.global_error("Error: all %d VCF files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d VCF files were missing or empty." % bad_file_count, continue_possible=True)

    #==========================================================================
    # Read in all vcf files and process into dict of SNPs passing various
    # criteria. Do this for each sample. Write to file.
    #==========================================================================
    if args.forceFlag or utils.target_needs_rebuild(list_of_vcf_files, snp_list_file_path):
        snp_dict = dict()
        excluded_sample_directories = set()
        for sample_dir, vcf_file_path in zip(sorted_list_of_sample_directories, list_of_vcf_files):

            if not os.path.isfile(vcf_file_path):
                continue
            if os.path.getsize(vcf_file_path) == 0:
                continue

            utils.verbose_print("Processing VCF file %s" % vcf_file_path)
            sample_name = os.path.basename(os.path.dirname(vcf_file_path))
            snp_set = utils.convert_vcf_file_to_snp_set(vcf_file_path)
            max_snps = args.maxSnps
            if max_snps >= 0 and len(snp_set) > max_snps:
                utils.verbose_print("Excluding sample %s having %d snps." % (sample_name, len(snp_set)))
                excluded_sample_directories.add(sample_dir)
                continue

            for key in snp_set:
                if key not in snp_dict:
                    sample_list = [sample_name]
                    snp_dict[key] = sample_list
                else:
                    sample_list = snp_dict[key]
                    sample_list.append(sample_name)

        utils.verbose_print('Found %d snp positions across %d sample vcf files.' % (len(snp_dict), len(list_of_vcf_files)))
        utils.write_list_of_snps(snp_list_file_path, snp_dict)
        utils.verbose_print("")

        #==========================================================================
        # Write the filtered list of sample directories
        #==========================================================================
        #sample_directories_list_path = sample_directories_list_path + ".filtered"
        filtered_sample_directories_list_path = args.filteredSampleDirsFile
        with open(filtered_sample_directories_list_path, "w") as filtered_samples_file_object:
            # Loop over the unsorted list to keep the order of samples the same as the original.
            # This will keep the same HPC log file suffix number.
            for sample_dir in unsorted_list_of_sample_directories:
                if sample_dir not in excluded_sample_directories:
                    filtered_samples_file_object.write("%s\n" % sample_dir)
    else:
        utils.verbose_print("SNP list %s has already been freshly built.  Use the -f option to force a rebuild." % snp_list_file_path)
    utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))


def call_consensus(args):
    """Call the consensus base for a sample

    Call the consensus base for a sample at the positions where SNPs were found
    in any of the samples.
    This function expects, or creates '(*)', the following
        files arranged in the following way:
            snplist.txt
            samples
                sample_name_one/reads.all.pileup
                sample_name_one/consensus.fasta (*)

    The files are used as follows:
        1. The snplist.txt input file contains the list of SNP positions
           extracted from all the var.flt.vcf files combined.
        2. The reads.all.pileup input file is a pileups at all positions
           used to determine the nucleotide base at each SNP position.
        3. The consensus.fasta output file contains the SNP calls for each
           sequence, arranged as a fasta file with one sequence per sample.

    The snplist.txt, and reads.all.pileup are created outside of this function.
       The package documentation provides an example
        of creating these files based on the lambda_virus sequence that is used
        as one test for this package.

    Args:
        forceFlag : boolean
            flag to force processing even when result file already exists and
            is newer than inputs
        snpListFile : str
            File path (not just file name) of text format list of SNP positions
        excludeFile : str
            File path of VCF file of positions to exclude from the snp matrix.
        allPileupFile : str
            Relative or absolute path to the genome-wide pileup file for this
            sample
        consensusFile : str
            Output file. Relative or absolute path to the consensus fasta file
            for this sample.
        minBaseQual : int
            Mimimum base quality score to count a read. All other snp filters
            take effect after the low-quality reads are discarded.
        minConsFreq : float
            Consensus frequency. Mimimum fraction of high-quality reads
            supporting the consensus to make a call.
        minConsStrdDpth : int
            Consensus strand depth. Minimum number of high-quality reads
            supporting the consensus which must be present on both the
            forward and reverse strands to make a call.
        minConsStrdBias : float
            Strand bias. Minimum fraction of the high-quality
            consensus-supporting reads which must be present on both the
            forward and reverse strands to make a call. The numerator of this
            fraction is the number of high-quality consensus-supporting reads
            on one strand.  The denominator is the total number of high-quality
            consensus-supporting reads on both strands combined.

    Raises:

    Examples:
    args = argparse.Namespace
    args.snpListFile = 'snplist.txt'
    args.allPileupFile = 'reads.all.pileup'
    args.consensusFile = 'consensus.fasta'
    args.minBaseQual = 15
    args.minConsFreq = 0.6
    args.minConsStrdDpth = 4
    args.minConsStrdBias = 0.10
    args.vcfFailedSnpGt = '.'
    call_consensus(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    snp_list_file_path = args.snpListFile
    all_pileup_file_path = args.allPileupFile
    sample_directory = os.path.dirname(os.path.abspath(all_pileup_file_path))
    sample_name = os.path.basename(sample_directory)
    consensus_file_path = args.consensusFile
    consensus_file_dir = os.path.dirname(os.path.abspath(consensus_file_path))
    vcf_file_name = args.vcfFileName
    vcf_file_path = os.path.join(consensus_file_dir, vcf_file_name) if vcf_file_name else None

    bad_file_count = utils.verify_existing_input_files("Snplist file", [snp_list_file_path])
    if bad_file_count > 0:
        utils.global_error("Error: cannot call consensus without the snplist file.")

    bad_file_count = utils.verify_non_empty_input_files("Pileup file", [all_pileup_file_path])
    if bad_file_count > 0:
        utils.sample_error("Error: cannot call consensus without the pileup file.", continue_possible=False)

    source_files = [snp_list_file_path, all_pileup_file_path]

    exclude_file_path = args.excludeFile
    if exclude_file_path:
        bad_file_count = utils.verify_existing_input_files("Exclude file", [exclude_file_path])
        if bad_file_count > 0:
            utils.sample_error("Error: cannot call consensus without the file of excluded positions.", continue_possible=False)
        excluded_positions = utils.convert_vcf_file_to_snp_set(exclude_file_path)
        source_files.append(exclude_file_path)
    else:
        excluded_positions = set()

    # Check if the result is already fresh
    if not args.forceFlag and not utils.target_needs_rebuild(source_files, consensus_file_path):
        utils.verbose_print("Consensus call file %s has already been freshly built.  Use the -f option to force a rebuild." % consensus_file_path)
        utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))
        return

    # Load the list of which positions to called
    snp_list = utils.read_snp_position_list(snp_list_file_path)
    snplist_length = len(snp_list)
    utils.verbose_print("snp position list length = %d" % snplist_length)
    utils.verbose_print("excluded snps list length = %d" % len(excluded_positions))
    utils.verbose_print("total snp position list length = %d" % (snplist_length + len(excluded_positions)))

    # Call consensus. Write results to file.
    position_consensus_base_dict = dict()

    caller = pileup.ConsensusCaller(args.minConsFreq,
                                    args.minConsStrdDpth,
                                    args.minConsStrdBias)

    snp_positions = set(snp_list)
    if args.vcfAllPos:
        parse_positions = None
    else:
        parse_positions = snp_positions.union(excluded_positions)
    pileup_reader = pileup.Reader(all_pileup_file_path,
                                  args.minBaseQual,
                                  parse_positions)
    if vcf_file_name:
        writer = vcf_writer.SingleSampleWriter(vcf_file_path, args.vcfPreserveRefCase)
        filters = caller.get_filter_descriptions()
        # TODO: it would be better if the exclude file contained filter headers we could read and re-use here instead of hard-coding this
        filters.append(("Region", "Position is in dense region of snps or near the end of the contig."))
        writer.write_header(sample_name, filters, args.vcfRefName)
    for pileup_record in pileup_reader:
        chrom = pileup_record.chrom
        pos = pileup_record.position
        consensus_base, fail_reasons = caller.call_consensus(pileup_record)
        if (chrom, pos) in excluded_positions:
            # TODO: it would be better if the exclude file contained filter reasons we could re-use here instead of hard coding this
            fail_reasons = fail_reasons or []
            fail_reasons.append("Region")
        if (chrom, pos) in snp_positions:
            if fail_reasons:
                position_consensus_base_dict[(chrom, pos)] = '-'
            else:
                position_consensus_base_dict[(chrom, pos)] = consensus_base

        if vcf_file_name:
            writer.write_from_pileup(pileup_record, fail_reasons, args.vcfFailedSnpGt)
    if vcf_file_name:
        writer.close()

    utils.verbose_print("called consensus positions = %i" % (len(position_consensus_base_dict)))

    consensus_list = [position_consensus_base_dict.get(key, '-') for key in snp_list]
    consensus_str = ''.join(consensus_list)
    snp_seq_record = SeqRecord(Seq(consensus_str), id=sample_name, description="")

    # Write the consensus calls to a fasta file
    with open(consensus_file_path, "w") as fasta_file_object:
        SeqIO.write([snp_seq_record], fasta_file_object, "fasta")

    utils.verbose_print("")
    utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))


def create_snp_matrix(args):
    """Create SNP matrix

    Description:
    Create the SNP matrix containing the consensus base for each of the samples
    at the positions where SNPs were found in any of the samples.  The matrix
    contains one row per sample and one column per SNP position.  Non-SNP
    positions are not included in the matrix.
    This function expects, or creates '(*)', the following
        files arranged in the following way:
            sampleDirectories.txt
            samples
                sample_name_one/consensus.fasta
                ...
            snpma.fasta (*)

    The files are used as follows:
        1. The sampleDirectories.txt input file contains a list of the paths to
           the sample directories.
        2. The consensus.fasta input files are previously called consensus
           for each sample to construct the SNP matrix fasta file.
        3. The snpma.fasta output file contains the SNP calls for each
           sequence, arranged as a multi-fasta file with one sequence per
           sample.

    The sampleDirectories.txt, and consensus.fasta are created outside of this
        function. The package documentation provides an example of creating
        these files based on the lambda_virus sequence that is used as one
        test for this package.

    Args:
        sampleDirsFile : str
            File path (not just file name) of file containing paths
            to directories containing consensus.fasta file for each sequence.
        snpListFile : str
            File path (not just file name) of text format list of SNP positions
        consFileName : str
            File name of the previously called consensus fasta files which must
            exist in each of the sample directories
        snpmaFile : str
            File path (not just file name) of the output snp matrix, formatted
            as a fasta file, with each sequence (all of identical length)
            corresponding to the SNPs in the correspondingly named sequence.

    Raises:

    Examples:
    args = argparse.Namespace
    args.sampleDirsFile = 'sampleDirectories.txt'
    args.consFileName = 'consensus.fasta'
    args.snpmaFile = 'snpma.fasta'
    args.minConsFreq = 0.6
    create_snp_matrix(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Prep work
    #==========================================================================
    sample_directories_list_filename = args.sampleDirsFile
    bad_file_count = utils.verify_non_empty_input_files("File of sample directories", [sample_directories_list_filename])
    if bad_file_count > 0:
        utils.global_error(None)

    with open(sample_directories_list_filename, "r") as sample_directories_list_file:
        list_of_sample_directories = [line.rstrip() for line in sample_directories_list_file]
    list_of_sample_directories = sorted([d for d in list_of_sample_directories if d])

    #==========================================================================
    # Verify input consensus.fasta files exist
    #==========================================================================
    consensus_files = []
    bad_file_count = 0
    for sample_directory in list_of_sample_directories:
        consensus_file_path = os.path.join(sample_directory, args.consFileName)
        bad_count = utils.verify_non_empty_input_files("Consensus fasta file", [consensus_file_path])
        if bad_count == 1:
            bad_file_count += 1
        else:
            consensus_files.append(consensus_file_path)  # keep the list of good files

    if bad_file_count == len(list_of_sample_directories):
        utils.global_error("Error: all %d consensus fasta files were missing or empty." % bad_file_count)
    elif bad_file_count > 0:
        utils.sample_error("Error: %d consensus fasta files were missing or empty." % bad_file_count, continue_possible=True)

    #==========================================================================
    # Check if the result is already fresh
    #==========================================================================
    snpma_file_path = args.snpmaFile
    source_files = consensus_files
    if not args.forceFlag:
        if not utils.target_needs_rebuild(source_files, snpma_file_path):
            utils.verbose_print("SNP matrix %s has already been freshly built.  Use the -f option to force a rebuild." % snpma_file_path)
            utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))
            return

    #==========================================================================
    #   Create snp matrix. Write results to file.
    #==========================================================================
    with open(snpma_file_path, "w") as output_file:
        for consensus_file_path in consensus_files:
            utils.verbose_print("Merging " + consensus_file_path)
            with open(consensus_file_path, "r") as input_file:
                for line in input_file:
                    output_file.write(line)

    utils.verbose_print("")
    utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))


def create_snp_reference_seq(args):
    """Write reference sequence bases at SNP locations to a fasta file.

    Description:
    Write reference sequence bases at SNP locations to a fasta file.
    This function expects, or creates '(*)', the following files:
            reference.fasta
            snplist.txt
            referenceSNP.fasta (*)

    The files are used as follows:
        1. The reference.fasta input file contains the whole-genome reference
           bases.
        2. The snplist.txt input file contains the list of SNP positions across
           all the samples.
        2. The referenceSNP.fasta output file contains the reference bases at
           the identified SNP locations.

    The snplist.txt file is created outside of this function.  The package
        documentation provides an example of creating this file based on the
        lambda_virus sequence that is used as one test for this package.

    Args:
        referenceFile: File path (not just file name) for reference sequence
            (in fasta format
        snpListFile: File path (not just file name) of text format list of SNP
            positions
        snpRefFile: File path (not just file name) for the SNP reference
            sequence file.

    Raises:

    Examples:
    args = argparse.Namespace
    args.referenceFile = 'reference.fasta'
    args.snpListFile = 'snplist.txt'
    args.snpRefFile = 'referenceSNP.fasta'
    create_snp_reference_seq(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    #    Write reference sequence bases at SNP locations to a fasta file.
    #==========================================================================
    reference_file = args.referenceFile
    snp_list_file_path = args.snpListFile
    snp_ref_seq_path = args.snpRefFile

    #==========================================================================
    # Verify input files exist
    #==========================================================================
    bad_file_count = utils.verify_existing_input_files("Snplist file", [snp_list_file_path])
    if bad_file_count > 0:
        utils.global_error("Error: cannot create the snp reference sequence without the snplist file.")

    bad_file_count = utils.verify_non_empty_input_files("Reference file", [reference_file])
    if bad_file_count > 0:
        utils.global_error("Error: cannot create the snp reference sequence without the reference fasta file.")

    #==========================================================================
    # Find the reference bases at the snp positions
    #==========================================================================
    source_files = [reference_file, snp_list_file_path]
    if args.forceFlag or utils.target_needs_rebuild(source_files, snp_ref_seq_path):
        utils.write_reference_snp_file(reference_file, snp_list_file_path, snp_ref_seq_path)
        utils.verbose_print("")
    else:
        utils.verbose_print("SNP reference sequence %s has already been freshly built.  Use the -f option to force a rebuild." % snp_ref_seq_path)

    utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))


def calculate_snp_distances(args):
    """Calculate pairwise sample SNP distances.

    Description:
    Calculate pairwise SNP distances from the multi-fasta SNP matrix.
    Generate a file of pairwise distances and a file containing a matrix
    of distances.
    This function expects, or creates '(*)', the following files:
            snpma.fasta
            snp_distance_pairwise.tsv*
            snp_distance_matrix.tsv*

    The files are used as follows:
        1. The snpma.fasta input file contains the snp matrix for all samples
        2. The snp_distance_pairwise.tsv output file contains a three column
            tab-separated table of distances between all pairs of samples
        2. The snp_distance_matrix.tsv output file contains a matrix of
           distances between all samples.

    Args:
        inputFile: File path (not just file name) for the snp matrix in fasta format
        pairwiseFile: File path (not just file name) of the output pairwise distance file
        matrixFile: File path (not just file name) for the output distance matrix file

    Raises:

    Examples:
    args = argparse.Namespace
    args.inputFile = 'snpma.fasta'
    args.pairwiseFile = 'snp_distance_pairwise.tsv'
    args.matrixFile = 'snp_distance_matrix.tsv'
    calculate_snp_distances(args)
    """
    utils.print_log_header()
    utils.print_arguments(args)

    #==========================================================================
    # Validate arguments
    #==========================================================================
    input_file = args.inputFile
    pairwise_file = args.pairwiseFile
    matrix_file = args.matrixFile
    force_flag = args.forceFlag

    bad_file_count = utils.verify_existing_input_files("SNP matrix file", [input_file])
    if bad_file_count > 0:
        utils.global_error("Error: cannot calculate sequence distances without the snp matrix file.")

    if not pairwise_file and not matrix_file:
        utils.global_error("Error: no output file specified.")

    #==========================================================================
    # Check freshness
    #==========================================================================
    rebuild_pairwise_file = pairwise_file and utils.target_needs_rebuild([input_file], pairwise_file)
    rebuild_matrix_file = matrix_file and utils.target_needs_rebuild([input_file], matrix_file)
    if force_flag or rebuild_pairwise_file or rebuild_matrix_file:

        #------------------------------
        # Read in snp matrix file
        #------------------------------
        seqs = {}
        with open(input_file) as ifile:
            for line in ifile:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    curr_sample = line.lstrip('>')
                    seqs[curr_sample] = ''
                else:
                    seqs[curr_sample] += str(line)

        #------------------------------
        # Count mismatches
        #------------------------------
        ids = sorted(seqs.keys())
        pairwise_mismatches = dict() # tuple (seq1 id, seq2 id) -> int

        for id1, id2 in itertools.combinations(ids, 2):
            mismatches = utils.calculate_sequence_distance(seqs[id1], seqs[id2])
            pairwise_mismatches[(id1, id2)] = mismatches
            pairwise_mismatches[(id2, id1)] = mismatches

        #------------------------------
        # Print distance files
        #------------------------------
        if pairwise_file:
            with open(pairwise_file, 'w') as p_out:
                p_out.write('%s\n' % '\t'.join(['Seq1', 'Seq2', 'Distance']))
                for id1, id2 in itertools.product(ids, ids):
                    mismatches = pairwise_mismatches.get((id1, id2), 0) # zero when id1=id2
                    p_out.write("%s\t%s\t%i\n" % (id1, id2, mismatches))

        if matrix_file:
            with open(matrix_file, 'w') as m_out:
                m_out.write('\t%s\n' % '\t'.join(ids)) # matrix header
                # write table of mismatches
                for id1 in ids:
                    mismatches = [pairwise_mismatches.get((id1, id2), 0) for id2 in ids]
                    mismatch_strs = map(str, mismatches)
                    m_out.write("%s\t%s\n" % (id1, '\t'.join(mismatch_strs)))

    else:
        utils.verbose_print("Distance files have already been freshly built.  Use the -f option to force a rebuild.")
    utils.verbose_print("# %s %s finished" % (utils.timestamp(), utils.program_name()))
