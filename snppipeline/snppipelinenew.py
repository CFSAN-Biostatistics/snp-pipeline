#!/usr/bin/env python2.7

from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import argparse
import os
import pprint
import utilsnew

def run_snp_pipeline(options_dict):
    """Create SNP matrix
     
    Description:
    Create a SNP matrix. This function expects or creates '(*)' the following
        files arranged in the following way:
            mainPath            
                reference.fasta
                path.txt
                snplist.txt (*)
                snpma.fasta (*)
            samples
                sample_name_one/reads.pileup (*)
                sample_name_one/var.flt.vcf
  
    The files are used as follows:
        1. The reference.fasta file is used for alignment of the sequence
            data in the running of the samtools pileup.
        2. The path.txt file contains a list of the paths to the sample
            directories.
        3. The snplist.txt file contains the list of SNPs extracted from the
            var.flt.vcf file.            
        4. The snpma.fasta file contains the SNP calls for each sequence,
            arranged as a fasta file with one sequence per sample.
        5. The reads.pileup files are used to determine the nucleotide base at
            each SNP position for each sample to construct the SNP fasta file.
        6. The variant file var.flt.vcf is used to construct the SNP position
            list. 
    
    The samtool pileups are run in parallel using the python multiprocessing
        package.
    
    The vcf files are created outside of this function. Here we provide an
        example of creating these files based on the lambda_virus sequence
        that we use as one test for this package:
        
        1. Align sequences to reference
            bowtie2-align -p 11 -q -x /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reference/lambda_virus -1 reads4_1.fq -2 reads4_2.fq > reads4.sam

        2. Convert to bam file with only mapped positions
            samtools view -bS -F 4 -o /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reads/reads1_F4.bam reads1.sam

        3. Convert to a sorted bam 
            samtools sort  reads4_F4.bam reads4_F4.sorted.bam

        4. To get a bcf file from the pileup and bam file
            samtools mpileup -uf /Users/james.pettengill/Downloads/bowtie2-2.2.0/example/reference/lambda_virus.fa reads1_F4.sorted.bam.bam | bcftools view -bvcg - > reads1_F4.bcf

        5. To convert bcf to vcf
            bcftools view reads1_F4.bcf | vcfutils.pl varFilter -D1000 > var1_F4.flt.vcf    

    Args:
        maxThread: (15) Max number of cocurrent threads.
        mainPath:  (no default) Directory containing all input files. Output
            files will also be written here.
        Reference: (no default) File name for reference sequence (in fasta
            format) for mapping.
        pathFileName: ("path.txt") Name of file containing full paths to
            directories containing information for each sequence.
        snplistFileName: Text format list of SNP positions in samples.
        snpmaFileName: File name for snp matrix, formatted as a fasta file,
            with each sequence (all of identical length) corresponding to the
            SNPs in the correspondingly named sequence.
        bamFileName: #TODO - do we actually ever use this?
        pileupFileName: Name for pileup files. One is generated for each
            sample, and placed in the corresponding directory for each sample.
     
    Raises:
 
    Examples:
    args_dict = {'maxThread':2,
                 'mainPath':'',
                 'Reference':'',
                 'pathFileName':'path.txt',
                 'snplistFileName':'snplist.txt',
                 'snpmaFileName':'snpma.fa',
                 'bamFileName':'reads.bam',
                 'pileupFileName':'reads.pileup'
                }
    run_snp_pipeline(options_dict) 
    """

    #==========================================================================
    #Prep work     
    #==========================================================================
    verbose = True
    verbose_print  = print         if verbose else lambda *a, **k: None
    verbose_pprint = pprint.pprint if verbose else lambda *a, **k: None

    sample_directories_list_filename = (options_dict['mainPath'] +
                                        options_dict['pathFileName'])
    list_of_sample_directories = [line.rstrip() for line in open(sample_directories_list_filename, "r")]
    #remove any blank rows that were read in
    list_of_sample_directories = filter(None,list_of_sample_directories)

    #==========================================================================
    #read in all vcf files and process into list of SNPs passing various
    #  criteria. Do this for each sample. Write to file
    #==========================================================================

    snp_list_dict = dict()
    
    for sample_directory in list_of_sample_directories:
        sample_name  = sample_directory.split(os.sep)[-1]
    
        #TODO - look at use of PyVCF to process vcf file
        for line in open(sample_directory+ "/var.flt.vcf","r"):
            curVcfFileLine=line.strip()
            if curVcfFileLine.startswith("#"):
                continue
            curLineData = curVcfFileLine.split()
            chrom = curLineData[0]
            pos   = curLineData[1]
            info  = curLineData[7]
            if str("INDEL") in info:
                continue
            infoFields = info.split(";")
            dpFlag = False
            af1Flag = False
            for infoField in infoFields:
                infoPair = infoField.split("=")
                if infoPair[0] == "DP" and int(infoPair[1]) >= 10:
                    dpFlag = True
                elif (infoPair[0] == "AF1" and infoPair[1] == "1" ) or (infoPair[0] == "AR" and infoPair[1] == "1.00"):
                    af1Flag = True
            # find a good record fo SNP position, save data to hash
            if dpFlag and af1Flag:
                if not snp_list_dict.has_key(chrom + "\t" + pos):
                    record = [1]
                    record.append(sample_name)
                    snp_list_dict[chrom + "\t" + pos] = record
                else:
                    record = snp_list_dict[chrom + "\t" + pos]
                    record[0] += 1
                    record.append(sample_name)
        
    snp_list_file_path=options_dict['mainPath'] + options_dict['snplistFileName']
    utilsnew.write_list_of_snps(snp_list_file_path,snp_list_dict)   
    
    #==========================================================================
    # Generate Pileups of samples (in parallel)
    # Note that we use map and not map_async so that we block
    #   until all the pileups are done (or bad things will happen in subsequent
    #   parts of the code).
    #==========================================================================
    
    #create a list of tuples containing values need for pileup code (as passed
    #  via pileup code wrapper)
    parameter_list = zip(list_of_sample_directories,
                         len(list_of_sample_directories)*[options_dict])
    
    verbose_print("Starting Pileups.")
    pool        = Pool(processes=options_dict['maxThread']) # start pool
    result_many = pool.map(utilsnew.pileup_wrapper, parameter_list) #parallel
    
    verbose_pprint(result_many)
    verbose_print("Pileups are finished.")
    
    #==========================================================================
    #   Create snp matrix. Write results to file.
    #==========================================================================
    
    records = []

    for sample_directory in list_of_sample_directories:
        sample_name       = sample_directory.split(os.sep)[-1]
        pileup_file_name  = sample_directory + "/reads.pileup"
        positionValueHash = utilsnew.create_consensus_dict(pileup_file_name)

        seqString = ""
        for key in sorted(snp_list_dict.iterkeys()):  #TODO - Why is sorting what we want to do?
            chrom,pos   = key.split()
            if positionValueHash.has_key(chrom + ":" + pos):
                seqString += positionValueHash[chrom + ":" + pos]
            else:
                seqString += "-"

        seq = Seq(seqString)
        seqRecord = SeqRecord(seq,id=sample_name)
        records.append(seqRecord)
    
    #Write bases for snps for each sequence to a fasta file           
    fastaFile = open(options_dict['mainPath'] + options_dict['snpmaFileName'], "w") 
    SeqIO.write(records, fastaFile, "fasta")
    fastaFile.close()
    
    #Write reference sequence bases at SNP locations to a fasta file
    if options_dict['includeReference']:
        snp_list_file_path       = options_dict['mainPath'] + options_dict['snplistFileName']
        reference_file_path      = options_dict['mainPath'] + options_dict['Reference']
        snp_reference_file_path  = options_dict['mainPath'] + "referenceSNP.fasta"   #TODO - should make this configurable
        utilsnew.write_reference_snp_file(reference_file_path,snp_list_file_path,snp_reference_file_path)


#==============================================================================
# Command line driver
#==============================================================================
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run SNP pipeline.')
    parser.add_argument('-n','--n-processes',      dest='maxThread',       type=int, default=4,help='Max number of concurrent jobs.')
    parser.add_argument('-d','--mainPath',         dest='mainPath',        type=str, default='/home/yan.luo/Desktop/analysis/Montevideo/XL-C2/bowtie/Matrices/',help='Path for all files')
    parser.add_argument('-r','--Reference',        dest='Reference',       type=str, default='reference.fasta',help='reference for mapping')
    parser.add_argument('-f','--pathFileName',     dest='pathFileName',    type=str, default='path.txt',help='Path file name')
    parser.add_argument('-l','--snplistFileName',  dest='snplistFileName', type=str, default='snplist.txt',help='Snplist file name')
    parser.add_argument('-a','--snpmaFileName',    dest='snpmaFileName',   type=str, default='snpma.fa',help='fasta file name')
    parser.add_argument('-b','--bamFileName',      dest='bamFileName',     type=str, default='reads.bam',help='bam file name')
    parser.add_argument('-p','--pileupFileName',   dest='pileupFileName',  type=str, default='reads.pileup',help='pileup file name')
    parser.add_argument('-v','--verbose',          dest='verbose',         type=int, default=1,help='Verbose flag (0=no info, 5=lots')
    parser.add_argument('-i','--includeReference', dest='includeReference',type=bool,default=False,help='Write reference sequence bases at SNP positions in fasta format.')
    parser.add_argument('-o','--useOldPileups',    dest='useOldPileups',   type=bool,default=False,help='Use available pileup files.')
    args_dict = vars(parser.parse_args())

    print("Running SNP pipeline with arguments:")
    pprint.pprint(args_dict)
    run_snp_pipeline(args_dict)


#----VCF file info to work on soon
#The header line names the 8 fixed, mandatory columns. These columns are as follows:
#
#    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,extra = string.split(current_line_data,maxsplit=9)
#
#If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number of sample IDs. The header line is tab-delimited.
#data_line = 'gi|9626243|ref|NC_001416.1|	43852	.	ACCT	A	214	.	INDEL;DP=40;VDB=0.0160;AF1=1;AC1=2;DP4=0,0,15,16;MQ=42;FQ=-128	GT:PL:GQ	1/1:255,93,0:99'
#gi|9626243|ref|NC_001416.1|	44530	.	GCAACA	GCA	214	.	INDEL;DP=48;VDB=0.0193;AF1=1;AC1=2;DP4=0,0,16,15;MQ=42;FQ=-128	GT:PL:GQ	1/1:255,93,0:99
#gi|9626243|ref|NC_001416.1|	46842	.	G	C	207	.	DP=41;VDB=0.0147;AF1=1;AC1=2;DP4=0,0,14,8;MQ=42;FQ=-90	GT:PL:GQ	1/1:240,63,0:99
#----

