"""
Utility functions for fastq files.
"""

from __future__ import print_function
import collections
import gzip
import os
import re


def sample_id(fastq_path):
    """Return the sample id of a fastq file as defined by the cfsan snp pipeline.

    This assumes the file has been placed in a directory whose basename is the sample_id.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file.

    Returns
    -------
    sample_id : str
        The cfsan snp pipeline defines the sample id as the basename of the parent directory.

    Examples
    --------
    >>> sample_id("x.fastq") == os.path.basename(os.getcwd())
    True

    >>> sample_id("aaa/x.fastq")
    'aaa'

    >>> sample_id("aaa/bbb/x.fastq")
    'bbb'
    """
    sample_abs_path = os.path.abspath(fastq_path)
    sample_dir = os.path.dirname(sample_abs_path)
    sample_id = os.path.basename(sample_dir)
    return sample_id


# Regular expression used to parse Illumina fastq sequence id lines.
# This is an abbreviated regular expression, there may be more tags after these.
ILLUMINA_FASTQ_SEQ_ID_REGEX = re.compile("@[^:]*(M[0-9]+):([^:]*):([^:]*):([^:]*)")

# Named tuple to contain fastq metadata
FastqSeqTags = collections.namedtuple("FastqSeqTags", "platform instrument flow_cell lane")

def parse_seqid_line(seqid_line):
    """Examine a fastq sequence id line and extract various metadata tags.

    Looks at a sequence id line of a fastq file to extract the instrument,
    flowcell, and lane if possible.  There is one of these lines preceeding
    every read in a fastq file and the information is highly redundant.

    Illumina fastq files have this sequence id format:
    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>

    Parameters
    ----------
    seqid_line : str
        A sequence identifier line from a fastq file.

    Returns
    -------
    FastqSeqTags named-tuple with the following named elements or None if the line cannot be parsed
        platform : str
            Currently, only Illumina data is supported and this will always be "illumina" if the line can be parsed
        instrument : str
            Instrument name
        flow_cell : str
            Flowcell with leading zeros removed
        lane : str
            Lane number

    Examples
    --------
    # Not Illumina format
    >>> parse_seqid_line("@SRR1206159_1/1") is None
    True

    # fastq-dump default settings
    >>> parse_seqid_line("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A1474', lane='1')

    # trailing colon after lane
    >>> parse_seqid_line("@M00229:7:A0WG8:1:")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A0WG8', lane='1')

    # no trailing colon after lane
    >>> parse_seqid_line("@M00229:7:A0WG8:1")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A0WG8', lane='1')

    # leading whitespace
    >>> parse_seqid_line("  @M00229:7:A0WG8:1")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A0WG8', lane='1')

    # hyphen before instrument and flowcell
    >>> parse_seqid_line("@-M00229:7:-A0WG8:1:")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A0WG8', lane='1')

    # Lee's fastq-dump defline-seq format
    >>> parse_seqid_line("@SRR498423_HWI-M00229:7:000000000-A0WG8:1:1:12203:2225/1")
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A0WG8', lane='1')

    References
    ----------
    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    """

    # Does the sequence id line look like Illumina format?
    match = ILLUMINA_FASTQ_SEQ_ID_REGEX.search(seqid_line)
    if match is None:
        return None

    # Parse the sequence id line
    instrument = match.group(1)
    flow_cell = match.group(3)
    lane = match.group(4)

    # strip leading zeros and minus from flowcell if present
    flow_cell_parts = flow_cell.split('-')
    flow_cell = flow_cell_parts[-1]

    return FastqSeqTags("illumina", instrument, flow_cell, lane)


def extract_metadata_tags(fastq_path):
    """Examine a fastq file and extract various metadata tags.

    Looks at the first line of the fastq file to extract the instrument,
    flowcell, and lane if possible.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file.

    Returns
    -------
    FastqSeqTags named-tuple with the following named elements or None if the line cannot be parsed
        platform : str
            Currently, only Illumina data is supported and this will always be "illumina" if the line can be parsed
        instrument : str
            Instrument name
        flow_cell : str
            Flowcell with leading zeros removed
        lane : str
            Lane number

    Examples
    --------
    # Setup tests
    >>> from tempfile import NamedTemporaryFile

    # Uncompressed file
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151\\n")
    >>> f.close()
    >>> extract_metadata_tags(filepath)
    FastqSeqTags(platform='illumina', instrument='M00229', flow_cell='A1474', lane='1')

    # Compressed file
    >>> os.rename(filepath, filepath + ".gz")
    >>> filepath = filepath + ".gz"
    >>> gf = gzip.open(filepath, "wt")
    >>> num_bytes = gf.write("@SRR498276.1 HWI-M00339:9:000000000-A1444:2:1:15012:1874 length=151\\n")
    >>> gf.close()
    >>> extract_metadata_tags(filepath)
    FastqSeqTags(platform='illumina', instrument='M00339', flow_cell='A1444', lane='2')

    >>> os.unlink(filepath)

    References
    ----------
    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    """
    # Case 1: file is gzipped
    if fastq_path.endswith(".gz"):
        line = ""
        f = gzip.open(fastq_path, "rt")
        try:
            line = f.readline()
        finally:
            f.close()
        return parse_seqid_line(line)

    # Case 2: file is NOT gzipped
    with open(fastq_path, "r") as f:
        line = ""
        line = f.readline()
    return parse_seqid_line(line)


# Named tuple to contain read group tags parsed from fastq metadata
ReadGroupTags = collections.namedtuple("ReadGroupTags", "ID SM LB PL PU")

def construct_read_group_tags(fastq_path):
    """Examine a fastq file and construct read group tags from the
    flowcell and lane if possible.

    A read group is a set of reads originating from a separate library
    and generated from a single run of a sequencing instrument.  This function
    assumes all the reads in the fastq file are from the same library.

    This function only looks at the first line of the fastq file and assumes
    all the reads in the file have the same read group.  When sequencing is
    spread over multiple lanes (NextSeq), the reads should be kept in separate
    fastq files until after the read groups are identified.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file.

    Returns
    -------
    tags : ReadGroupTags or None
        ID : str
            An identifier for the source of reads within the fastq file formed
            from the flowcell and lane.  The read group id is {flowcell}.{lane}".
        SM : str
            The name of the sample sequenced in this read group.
        LB : str
            DNA preparation library identifier.  This function assumes all the reads
            in the fastq file are from the same library.  Always set to "1".
        PL : str
            Platform/technology used to produce the read.  Valid values are:
            illumina, solid, ls454, helicos and pacbio.
        PU : str
            Platform Unit is {flowcell}.{lane}.{sample}.

    Examples
    --------
    # Setup tests
    >>> from tempfile import NamedTemporaryFile

    # Sequence id line contains flowcell and lane
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 HWI-M00229:9:000000000-A1474:1:1:15012:1874 length=151\\n")
    >>> f.close()
    >>> construct_read_group_tags(filepath)
    ReadGroupTags(ID='A1474.1', SM='tmp', LB='1', PL='illumina', PU='A1474.1.tmp')
    >>> os.unlink(filepath)

    # Sequence id line missing flowcell and lane
    >>> f = NamedTemporaryFile(delete=False, mode='w')
    >>> filepath = f.name
    >>> num_bytes = f.write("@SRR498276.1 length=151\\n")
    >>> f.close()
    >>> construct_read_group_tags(filepath) is None
    True
    >>> os.unlink(filepath)

    References
    ----------
    http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    """
    tags = extract_metadata_tags(fastq_path)
    if tags is None:
        return None

    id = tags.flow_cell + '.' + tags.lane
    sm = sample_id(fastq_path)
    lb = '1'
    pl = tags.platform
    pu = id + '.' + sm
    return ReadGroupTags(id, sm, lb, pl, pu)
