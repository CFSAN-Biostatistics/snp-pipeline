"""
VCF file writer.  Typical workflow would follow a pattern like this:

import vcf_writer
import pileup

caller = pileup.ConsensusCaller(min_freq, min_strand_depth, min_strand_bias)
pileup_reader = pileup.Reader(file, min_base_quality, chrom_position_set=None)
writer = vcf_writer.SingleSampleWriter(out_file)
filters = caller.get_filter_descriptions()
writer.write_header(sample_id, filters, reference)
for pileup_record in pileup_reader:
    alt, fail_reasons = caller.call_consensus(record)
    writer.write_from_pileup(pileup_record, fail_reasons)
writer.close()
"""

from __future__ import print_function
from __future__ import absolute_import
import collections
import datetime
import sys
if sys.version_info < (3,):
    from StringIO import StringIO
else:
    from io import StringIO
import vcf
from snppipeline.__init__ import __version__
from snppipeline import pileup

VCF_VERSION = '##fileformat=VCFv4.1\n'
VCF_DATE    = '##fileDate=%Y%m%d\n'
VCF_SOURCE  = '##source=CFSAN SNP-Pipeline %s\n' % __version__
VCF_INFO    = '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">\n'
VCF_FILTER  = '##FILTER=<ID=%s,Description="%s">\n'
VCF_FORMAT  = '''
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filters using the same codes as the FILTER data element">
'''
VCF_REFERENCE = '##reference=%s\n'
VCF_HDR_LINE  = '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  %s\n'


#==============================================================================
#Define classes and functions
#==============================================================================

class SingleSampleWriter(object):
    """
    Class to write VCF files for a single sample given records
    parsed by the pileup module.
    """

    def __init__(self, file_path):
        """
        Initialize the VCF writer.
        """
        self.file_path = file_path
        self.file_handle = open(file_path, "w")

    def close(self):
        """
        Close the VCF writer.
        This method should be called after iterating through all the records.
        """
        if self.pyvcf_writer:
            self.pyvcf_writer.close()
        self.file_handle.close()

    def write_header(self, sample_id, filters, reference):
        """
        Write the VCF file header with the standard SNP Pipeline data elements.

        Parameters
        ----------
        sample_id : str
            Sample ID which will be written to the header line.
        filters : list of tuple(str, str)
            List of names and descriptions of filters which will be combined 
            and written to the header filter lines.
        reference : str
            Reference name which will be written to the header reference line.
        """
        # Write the template header to an in-memory buffer
        in_memory_file = StringIO()
        in_memory_file.name = "header.vcf"
        in_memory_file.write(VCF_VERSION)
        in_memory_file.write(datetime.datetime.strftime(datetime.datetime.now(), VCF_DATE))
        in_memory_file.write(VCF_SOURCE)
        in_memory_file.write(VCF_INFO)
        in_memory_file.write(VCF_FILTER % ("PASS", "All filters passed"))
        for name, description in filters:
            in_memory_file.write(VCF_FILTER % (name, description))
        in_memory_file.write(VCF_FORMAT)
        in_memory_file.write(VCF_REFERENCE % reference)
        in_memory_file.write(VCF_HDR_LINE  % sample_id)

        # Rewind to the beginning of the file buffer to prepare for reading
        in_memory_file.seek(0)

        # Feed the template to pyVcf and write the header to our vcf file
        vcf_template = vcf.Reader(in_memory_file)
        self.pyvcf_writer = vcf.Writer(self.file_handle, template=vcf_template)

        # Extract the format string from the header.  It will be the same for
        # all positions, so only do this once.
        format_lines = VCF_FORMAT.split('\n')
        format_lines = [line for line in format_lines if len(line) > 0]
        format_lines = [line.replace("##FORMAT=<ID=", "") for line in format_lines]
        tokens = [line.split(',')[0] for line in format_lines]
        self.format_str = ':'.join(tokens)
        self.VcfCallData = collections.namedtuple('VcfCallData', tokens)  # this creates a new class called VcfCallData

    def _make_vcf_record_from_pileup(self, pileup_record, failed_filters):
        """
        Create a PyVCF model record for a single sample from a single pileup 
        record.

        Parameters
        ----------
        pileup_record : pileup.Record
            Previously parsed pileup record.
        failed_filters : list of str or None
            List of failed snp filters or None if all filters passed

        Returns
        -------
        record : vcf.model._Record
            PyVCF record with attributes describing chrom, pos, ref, alt, depth, etc.

        Examples
        --------
        >>> writer = SingleSampleWriter("/dev/null")
        >>> caller = pileup.ConsensusCaller(0.5, 0, 0.0) # freq, strand_depth, strand_bias
        >>> filters = caller.get_filter_descriptions()
        >>> writer.write_header("dummy-sample-name", filters, "dummy_reference-name")
        >>>
        >>> r = pileup.Record(['ID', 42, 'G', 14, 'aaaaAAAA...,,,', '00001111222333'], 15)
        >>> v = writer._make_vcf_record_from_pileup(r, None)
        >>> v.CHROM, v.POS, v.ID, v.REF, v.ALT, v.QUAL, v.FILTER, v.INFO, v.FORMAT
        ('ID', 42, None, 'G', 'A', None, [], {'NS': 1}, 'GT:SDP:RD:AD:RDF:RDR:ADF:ADR:FT')
        >>> v.samples[0]
        Call(sample=None, VcfCallData(GT='1', SDP=14, RD=6, AD=8, RDF=3, RDR=3, ADF=4, ADR=4, FT='PASS'))
        >>>
        >>> r = pileup.Record(['ID', 42, 'g', 14, 'aaaaAAAA...,,,', '00001111222333'], 15)
        >>> v = writer._make_vcf_record_from_pileup(r, None)
        >>> v.CHROM, v.POS, v.ID, v.REF, v.ALT, v.QUAL, v.FILTER, v.INFO, v.FORMAT
        ('ID', 42, None, 'g', 'A', None, [], {'NS': 1}, 'GT:SDP:RD:AD:RDF:RDR:ADF:ADR:FT')
        >>> v.samples[0]
        Call(sample=None, VcfCallData(GT='1', SDP=14, RD=6, AD=8, RDF=3, RDR=3, ADF=4, ADR=4, FT='PASS'))
        >>>
        >>> r = pileup.Record(['ID', 42, 'G', 14, 'gaaaGGGG...,,,', '00001111222333'], 15)
        >>> v = writer._make_vcf_record_from_pileup(r, None)
        >>> v.CHROM, v.POS, v.ID, v.REF, v.ALT, v.QUAL, v.FILTER, v.INFO, v.FORMAT
        ('ID', 42, None, 'G', '.', None, [], {'NS': 1}, 'GT:SDP:RD:AD:RDF:RDR:ADF:ADR:FT')
        >>> v.samples[0]
        Call(sample=None, VcfCallData(GT='0', SDP=14, RD=11, AD=3, RDF=7, RDR=4, ADF=0, ADR=3, FT='PASS'))
        >>>
        >>> r = pileup.Record(['ID', 42, 'g', 14, 'ggggGGGG...,,,', '00001111222333'], 15)
        >>> v = writer._make_vcf_record_from_pileup(r, None)
        >>> v.CHROM, v.POS, v.ID, v.REF, v.ALT, v.QUAL, v.FILTER, v.INFO, v.FORMAT
        ('ID', 42, None, 'g', '.', None, [], {'NS': 1}, 'GT:SDP:RD:AD:RDF:RDR:ADF:ADR:FT')
        >>> v.samples[0]
        Call(sample=None, VcfCallData(GT='0', SDP=14, RD=14, AD=0, RDF=7, RDR=7, ADF=0, ADR=0, FT='PASS'))
        >>> writer.close()
        """
        ref = pileup_record.reference_base
        upper_ref = ref.upper()
        alt = pileup_record.most_common_base

        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth">
        ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
        ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
        ##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
        ##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
        ##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
        ##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
        ##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filters using the same codes as the FILTER data element">

        if alt is None:  # no good depth 
            gt = '.'
            alt = '.'
            ad = 0
            adf = 0
            adr = 0
        elif alt == upper_ref:
            # When we do not call an alt, the alt depths are counted as anything not matching the reference
            gt = '0'
            alt = '.'
            ad = pileup_record.good_depth - pileup_record.base_good_depth[upper_ref]
            adf = pileup_record.forward_good_depth - pileup_record.forward_base_good_depth[upper_ref]
            adr = pileup_record.reverse_good_depth - pileup_record.reverse_base_good_depth[upper_ref]
        else:
            # When we call an alt, the alt depths are only the specific alt base called
            gt = '1'
            ad = pileup_record.base_good_depth[alt]
            adf = pileup_record.forward_base_good_depth[alt]
            adr = pileup_record.reverse_base_good_depth[alt]
        sdp = pileup_record.raw_depth
        rd = pileup_record.base_good_depth[upper_ref]
        rdf = pileup_record.forward_base_good_depth[upper_ref]
        rdr = pileup_record.reverse_base_good_depth[upper_ref]
        ft = ';'.join(failed_filters) if failed_filters else "PASS"

        #print("Position:", pileup_record.position, "Alt:", alt, gt, sdp, rd, ad, rdf, rdr, adf, adr, ft)
        sample_data = self.VcfCallData(gt, sdp, rd, ad, rdf, rdr, adf, adr, ft)

        # http://pyvcf.readthedocs.org/en/latest/API.html#vcf-model-call
        # vcf.model._Call(site, sample, data)
        call = vcf.model._Call(None, None, sample_data)

        # info_dict contains fields that match our VCF_HEAD info fields
        # this sample has a pileup record at this position, so NS=1
        info_dict = {'NS' : 1}

        # https://pyvcf.readthedocs.org/en/latest/API.html#vcf-model-record
        # vcf.model._Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes, samples=None)
        vcf_record = vcf.model._Record(
            pileup_record.chrom, 
            pileup_record.position, 
            None,  # Id
            ref,
            alt or '.',
            None, # Qual
            failed_filters or [],
            info_dict,
            self.format_str, # Format
            None, # Sample indexes
            [call]  # Samples
        )
        return vcf_record


    def write_from_pileup(self, pileup_record, failed_filters):
        """
        Write a single VCF record for a single sample from a single pileup 
        record.

        Parameters
        ----------
        pileup_record : pileup.Record
            Previously parsed pileup record.
        failed_filters : list of str or None
            List of failed snp filters or None if all filters passed
        """
        # Create the PyVCF record
        vcf_record = self._make_vcf_record_from_pileup(pileup_record, failed_filters)

        # Tell PyVCF to write the record
        self.pyvcf_writer.write_record(vcf_record)
