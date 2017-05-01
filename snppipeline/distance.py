"""This module is part of the CFSAN SNP Pipeline. It contains the code to
calculate the pairwise SNP distances between samples.
"""

from __future__ import print_function
from __future__ import absolute_import

import itertools

from snppipeline import utils
from snppipeline.utils import verbose_print


def calculate_snp_distances(args):
    """Calculate pairwise sample SNP distances.

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

    Parameters
    ----------
    args : Namespace
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
        verbose_print("# %s %s" % (utils.timestamp(), "Calculating all pairwise distances"))
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
