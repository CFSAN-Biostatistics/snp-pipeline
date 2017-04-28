"""This module is part of the CFSAN SNP Pipeline. It contains the code to
call the consensus base for a sample from a pileup file.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

from snppipeline import pileup
from snppipeline import utils
from snppipeline import vcf_writer


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

    Parameters
    ----------
    args : namespace
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
