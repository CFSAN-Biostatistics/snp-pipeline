"""This module is part of the CFSAN SNP Pipeline. It contains the code to
calculate quality metrics and SNP metrics for a single sample.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
import os
from tempfile import NamedTemporaryFile
import vcf

from snppipeline import command
from snppipeline import fastq
from snppipeline import utils
from snppipeline.utils import verbose_print


error_list = []


def handle_error(message):
    """Print an error message and store it in the list of errors.
    """
    global error_list

    verbose_print(message)
    error_list.append(message)


def verify_input_file(error_prefix, file_path):
    """Verify a file exists and is not empty.  Errors are logged and
    also captured in the metrics.tsv output.

    NOTE: there is a similar function in the utils module, but this is
    different and has special behavior just for metrics collection.  Here,
    the errors will not stop execution and logged differently.

    Parameters
    ----------
    error_prefix : str
        First part of error message to be logged if the file fails verifcation
    file_path : str
        Relative or absolute paths to file

    Returns
    -------
    ok : bool
        True if the file exists and is not empty
    """
    file_base_name = os.path.basename(file_path)
    if not os.path.isfile(file_path):
        handle_error(error_prefix + ' ' + file_base_name + " was not found.")
        return False
    if os.path.getsize(file_path) == 0:
        handle_error(error_prefix + ' ' + file_base_name + " is empty.")
        return False
    return True


def count_vcf_file_snps(file_path):
    """Scan a VCF file and count the number of snps.  If the VCF records
    have FT data elements, any value other than PASS will be ignored and
    not counted as a snp.

    Parameters
    ----------
    file_path : str
        Relative or absolute paths to file

    Returns
    -------
    num_snps : int or None
        Number of snps in the VCF file if there are any records.
        None if there are no records in the VCF file.
    """
    num_records = 0
    num_snps = 0
    with open(file_path) as f:
        reader = vcf.VCFReader(f)
        for record in reader:
            num_records += 1
            if not record.is_snp: # is ALT not in [A,C,G,T,N,*] ?
                continue
            if record.ALT == record.REF:
                continue
            for sample in record.samples:
                if not sample.is_variant: # is GT == REF ?
                    continue
                try:
                    if sample.data.FT != "PASS":
                        continue # don't count failed snps
                except AttributeError:
                    pass # if there is no FT element, we should assume the snp passed and count it
                num_snps += 1
    return num_snps


def count_missing_snp_matrix_positions(file_path, sample_id):
    """Count the number of missing positions in a fasta snp matrix file.

    Parameters
    ----------
    file_path : str
        Relative or absolute paths to file
    sample_id : str
        Only count the gaps for records with id matching this string

    Returns
    -------
    num_missing : int
        Number of gaps '-' in the fasta file
    """
    with open(file_path, "rU") as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == sample_id:
                return record.seq.count('-')
    return 0


def collect_metrics(args):
    """Collect the quality metrics and SNP metrics for a sample.

    This function expects, or creates '(*)', the following files arranged
    in the following way:
            reference
                referenceFile.fasta
            samples
                sample_name_one/*.fastq.gz
                sample_name_one/reads.sam
                sample_name_one/reads.sorted.deduped.bam
                sample_name_one/reads.sorted.bam
                sample_name_one/reads.all.pileup
                sample_name_one/var.flt.vcf
                sample_name_one/var.flt_preserved.vcf
                sample_name_one/consensus.fasta
                sample_name_one/consensus_preserved.fasta
                sample_name_one/consensus.vcf
                sample_name_one/consensus_preserved.vcf
                sample_name_one/metrics*

    The input files are created outside of this function. The package
    documentation provides an example of preparing these files based on the
    lambda_virus sequence that is used as one test for this package.

    Parameters
    ----------
    args : argparse.Namespace
        referenceFile : File path of the reference fasta file
        sampleDir : Relative or absolute directory of the sample
        consensusFastaFileName : File name of the consensus fasta file which must exist in the sample directory
        consensusPreservedFastaFileName : File name of the consensus preserved fasta file which must exist in the sample directory
        consensusVcfFileName : File name of the consensus vcf file which must exist in the sample directory
        consensusPreservedVcfFileName : File name of the consensus preserved vcf file which must exist in the sample directory
        maxSnps : Maximum allowed number of SNPs per sample
        metricsFile : Output file.  Relative or absolute path to the metrics file
    """
    utils.print_log_header(classpath=True)
    utils.print_arguments(args)

    #==========================================================================
    # Validate inputs
    #==========================================================================

    # Verify reference fasta file exists and is not empty
    reference_file_path = args.referenceFile
    utils.verify_non_empty_input_files("Reference file", [reference_file_path], error_handler="global")

    sample_dir = args.sampleDir
    utils.verify_non_empty_directory("Sample directory", sample_dir, error_handler="sample", continue_possible=False)

    metrics_file_path = args.metricsFile
    max_allowed_snps = args.maxSnps
    consensus_vcf_file_name = args.consensusVcfFileName
    consensus_preserved_vcf_file_name = args.consensusPreservedVcfFileName
    consensus_fasta_file_name = args.consensusFastaFileName
    consensus_preserved_fasta_file_name = args.consensusPreservedFastaFileName

    sample_id = utils.sample_id_from_dir(sample_dir)

    #==========================================================================
    # Read existing metrics file so some metrics can be reused
    #==========================================================================
    try:
        metrics = utils.read_properties(metrics_file_path)
    except IOError:
        metrics = dict()

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Get machine and flowcell from fastq header"))
    #-------------------------
    machine = ""
    flowcell = ""
    fastq_files = fastq.list_fastq_files(sample_dir)
    fastq_files = [f for f in fastq_files if os.path.isfile(f)] # Exclude broken symlinks
    if not fastq_files:
        handle_error("No fastq files were found.")
    else:
        tags = fastq.extract_metadata_tags(fastq_files[0])
        if tags:
            machine = tags.instrument or ""
            flowcell = tags.flow_cell or ""

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Sum file sizes of paired fastq files"))
    #-------------------------
    fastq_file_size = ""
    fastq_file_list = ""
    if fastq_files:
        fastq_file_size = sum([os.path.getsize(file) for file in fastq_files])

    # Make a comma separated list of just the fastq file names without directories
    fastq_file_list = [os.path.basename(file) for file in fastq_files]
    fastq_file_list = ", ".join(fastq_file_list)

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Calculate number of reads and %mapped from sam file"))
    #-------------------------
    num_reads = ""
    percent_reads_mapped = ""
    file = os.path.join(sample_dir, "reads.sam")
    if verify_input_file("SAM file", file):
        # Metrics already freshly collected?
        needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
        if not args.forceFlag and not needs_rebuild:
            num_reads = metrics.get("numberReads", "") # reuse already fresh metrics
            percent_reads_mapped = metrics.get("percentReadsMapped", "") # reuse already fresh metrics
        if num_reads and percent_reads_mapped:
            verbose_print("Reusing previously calculated number of reads and %mapped")
        else:
            num_reads = command.run("samtools view -S -c " + file)
            num_reads = num_reads.strip()
            mapped = command.run("samtools view -S -c -F 4 " + file)
            mapped = mapped.strip()
            try:
                percent_reads_mapped = 100.0 * float(mapped) / float(num_reads)
                percent_reads_mapped = "%.2f" % percent_reads_mapped
            except ValueError:
                handle_error("Cannot calculate number of reads and %mapped.")

    #-------------------------
    # Calculate number of duplicate reads from deduped bam file
    #-------------------------
    num_dup_reads = ""
    remove_duplicate_reads = os.environ.get("SnpPipeline_RemoveDuplicateReads") or "true"
    remove_duplicate_reads = remove_duplicate_reads.lower()
    if remove_duplicate_reads == "true":
        verbose_print("# %s %s" % (utils.timestamp(), "Calculate number of duplicate reads from deduped bam file"))
        file = os.path.join(sample_dir, "reads.sorted.deduped.bam")
        if verify_input_file("Deduped BAM file", file):
            # Metrics already freshly collected?
            needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
            if not args.forceFlag and not needs_rebuild:
                num_dup_reads = metrics.get("numberDupReads", "") # reuse already fresh metrics
            if num_dup_reads:
                verbose_print("Reusing previously calculated number of duplicate reads")
            else:
                num_dup_reads = command.run("samtools view -S -c -f 1024 " + file)
                num_dup_reads = num_dup_reads.strip()

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Calculate mean insert size from bam file"))
    #-------------------------
    ave_insert_size = ""
    file = os.path.join(sample_dir, "reads.sorted.bam")
    if verify_input_file("BAM file", file):
        # Metrics already freshly collected?
        needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
        if not args.forceFlag and not needs_rebuild:
            ave_insert_size = metrics.get("aveInsertSize", "") # reuse already fresh metrics
        if ave_insert_size:
            verbose_print("Reusing previously calculated mean insert size")
        else:
            # Extract inferred insert sizes (TLEN, column 9 of BAM file) for reads "mapped in proper pair" (2) and "first in pair" (64) = 66
            tempfile = NamedTemporaryFile(delete=False, dir=sample_dir, prefix="tmp.inserts.", mode='w')
            command.run("samtools view -f 66 " + file + " | cut -f 9 | sed 's/^-//'", tempfile.name)
            insert_count = 0
            insert_sum = 0
            with open(tempfile.name) as f:
                for line in f:
                    try:
                        insert_sum += int(line)
                        insert_count += 1
                    except ValueError:
                        pass
            os.unlink(tempfile.name)
            if insert_count > 0 and insert_sum > 0:
                ave_insert_size = float(insert_sum) / float(insert_count)
                ave_insert_size = "%.2f" % ave_insert_size
            else:
                handle_error("Cannot calculate mean insert size.")

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Calculate mean depth from pileup file"))
    #-------------------------
    ave_pileup_depth = ""
    file = os.path.join(sample_dir, "reads.all.pileup")
    if verify_input_file("Pileup file", file):
        # Metrics already freshly collected?
        needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
        if not args.forceFlag and not needs_rebuild:
            ave_pileup_depth = metrics.get("avePileupDepth", "") # reuse already fresh metrics
        if ave_pileup_depth:
            verbose_print("Reusing previously calculated mean pileup depth")
        else:
            depth_sum = 0
            with open(file) as f:
                for line in f:
                    tokens = line.split()
                    try:
                        depth_sum += int(tokens[3])
                    except (ValueError, IndexError):
                        pass
            reference_length = 0
            for record in SeqIO.parse(reference_file_path, "fasta"):
                reference_length += len(record)
            if depth_sum > 0 and reference_length > 0:
                #print("depth_sum=%i" % depth_sum);
                #print("reference_length=%i" % reference_length)
                ave_pileup_depth = float(depth_sum) / float(reference_length)
                ave_pileup_depth = "%.2f" % ave_pileup_depth
            else:
                handle_error("Cannot calculate mean pileup depth.")

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count number of high confidence SNP positions from phase 1 vcf file"))
    #-------------------------
    phase1_snps = ""
    excluded_sample = ""
    file = os.path.join(sample_dir, "var.flt.vcf")
    if verify_input_file("VCF file", file):
        # Metrics already freshly collected?
        needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
        if not args.forceFlag and not needs_rebuild:
            phase1_snps = metrics.get("phase1Snps", "") # reuse already fresh metrics
        if phase1_snps:
            verbose_print("Reusing previously calculated phase1 snps")
        else:
            phase1_snps = count_vcf_file_snps(file)

        # Flag excessive snps
        if max_allowed_snps > 0 and phase1_snps > max_allowed_snps:
            excluded_sample = "Excluded"
            handle_error("Excluded: exceeded %i maxsnps." % max_allowed_snps)
        phase1_snps = str(phase1_snps)

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count number of snp_filter preserved high confidence SNP positions from phase 1 vcf file"))
    #-------------------------
    phase1_snps_preserved = ""
    excluded_sample_preserved = ""
    file = os.path.join(sample_dir, "var.flt_preserved.vcf")
    if verify_input_file("VCF file", file):
        # Metrics already freshly collected?
        needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
        if not args.forceFlag and not needs_rebuild:
            phase1_snps_preserved = metrics.get("phase1SnpsPreserved", "") # reuse already fresh metrics
        if phase1_snps_preserved:
            verbose_print("Reusing previously calculated preserved phase1 snps")
        else:
            phase1_snps_preserved = count_vcf_file_snps(file)

        # Flag excessive snps
        if max_allowed_snps > 0 and phase1_snps_preserved > max_allowed_snps:
            excluded_sample_preserved = "Excluded"
            handle_error("Excluded: preserved exceeded %i maxsnps." % max_allowed_snps)
        phase1_snps_preserved = str(phase1_snps_preserved)

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count number of consensus snps from consensus vcf file"))
    #-------------------------
    phase2_snps = ""
    file = os.path.join(sample_dir, consensus_vcf_file_name)
    if verify_input_file("Consensus VCF file", file):
        # Omit the phase2 snp count if the sample is excluded.
        # It will be meaningless since this sample's phase1 snps are excluded from the snplist.
        if excluded_sample != "Excluded":
            # Metrics already freshly collected?
            needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
            if not args.forceFlag and not needs_rebuild:
                phase2_snps = metrics.get("snps", "") # reuse already fresh metrics
            if phase2_snps:
                verbose_print("Reusing previously calculated phase2 snps")
            else:
                phase2_snps = count_vcf_file_snps(file)
                phase2_snps = str(phase2_snps)

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count number of preserved consensus snps from consensus vcf file"))
    #-------------------------
    phase2_snps_preserved = ""
    file = os.path.join(sample_dir, consensus_preserved_vcf_file_name)
    if verify_input_file("Consensus VCF file", file):
        # Omit the phase2 snp count if the sample is excluded.
        # It will be meaningless since this sample's phase1 snps are excluded from the snplist.
        if excluded_sample_preserved != "Excluded":
            # Metrics already freshly collected?
            needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
            if not args.forceFlag and not needs_rebuild:
                phase2_snps_preserved = metrics.get("snpsPreserved", "") # reuse already fresh metrics
            if phase2_snps_preserved:
                verbose_print("Reusing previously calculated preserved phase2 snps")
            else:
                phase2_snps_preserved = count_vcf_file_snps(file)
                phase2_snps_preserved = str(phase2_snps_preserved)

    #------------------------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count missing positions in the snp matrix"))
    #------------------------------------------
    missing_pos = ""
    file = os.path.join(sample_dir, consensus_fasta_file_name)
    if verify_input_file("Consensus fasta file", file):
        # Omit the phase2 gap count if the sample is excluded.
        # It will be meaningless since this sample's phase1 snps are excluded from the snplist.
        if excluded_sample != "Excluded":
            # Metrics already freshly collected?
            needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
            if not args.forceFlag and not needs_rebuild:
                missing_pos = metrics.get("missingPos", "") # reuse already fresh metrics
            if missing_pos:
                verbose_print("Reusing previously calculated missing positions")
            else:
                missing_pos = count_missing_snp_matrix_positions(file, sample_id)
                missing_pos = str(missing_pos)

    #------------------------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Count missing positions in the preserved snp matrix"))
    #------------------------------------------
    missing_pos_preserved = ""
    file = os.path.join(sample_dir, consensus_preserved_fasta_file_name)
    if verify_input_file("Consensus fasta file", file):
        # Omit the phase2 gap count if the sample is excluded.
        # It will be meaningless since this sample's phase1 snps are excluded from the snplist.
        if excluded_sample_preserved != "Excluded":
            # Metrics already freshly collected?
            needs_rebuild = utils.target_needs_rebuild([file], metrics_file_path)
            if not args.forceFlag and not needs_rebuild:
                missing_pos_preserved = metrics.get("missingPosPreserved", "") # reuse already fresh metrics
            if missing_pos_preserved:
                verbose_print("Reusing previously calculated missing positions")
            else:
                missing_pos_preserved = count_missing_snp_matrix_positions(file, sample_id)
                missing_pos_preserved = str(missing_pos_preserved)

    #-------------------------
    verbose_print("# %s %s" % (utils.timestamp(), "Print results"))
    #-------------------------
    with open(metrics_file_path, "w") as f:
        print("sample=" + '"' + sample_id + '"', file=f)
        print("fastqFileList=" + '"' + fastq_file_list + '"', file=f)
        print("fastqFileSize=" + str(fastq_file_size), file=f)
        print("machine=" + machine, file=f)
        print("flowcell=" + flowcell, file=f)
        print("numberReads=" + num_reads, file=f)
        print("numberDupReads=" + num_dup_reads, file=f)
        print("percentReadsMapped=" + percent_reads_mapped, file=f)
        print("aveInsertSize=" + ave_insert_size, file=f)
        print("avePileupDepth=" + ave_pileup_depth, file=f)
        print("phase1Snps=" + phase1_snps, file=f)
        print("phase1SnpsPreserved=" + phase1_snps_preserved, file=f)
        print("snps=" + phase2_snps, file=f)
        print("snpsPreserved=" + phase2_snps_preserved, file=f)
        print("missingPos=" + missing_pos, file=f)
        print("missingPosPreserved=" + missing_pos_preserved, file=f)
        print("excludedSample=" + excluded_sample, file=f)
        print("excludedSamplePreserved=" + excluded_sample_preserved, file=f)
        print("errorList=" + '"' + ' '.join(error_list) + '"', file=f)
