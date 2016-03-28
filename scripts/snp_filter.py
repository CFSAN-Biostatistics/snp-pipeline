import sys
import os
import vcf
import argparse


#Define functions
"""
This function reads the SNP matrix from a SNP matrix file, filter SNPs with abnormal density (which are exported into individual VCF files, 
one VCF file per sample, and the file name is in format "SampleName_FilteredSNP.vcf"), and write other SNPs into a new SNP Matrix file,with 
the file name "OldFileName_Filtered". 
OldSNPMatrixFile: the file name of the original SNP matrix.

def FilterDenseSNPsFromMatrix(OldSNPMatrixFile):
    with open(OldSNPMatrixFile, 'r') as f1:
        line=f1.readline()
        if (len(line) == 0):
            #End of the old SNP matrix file
            ProcessSNPs()
            return
        sa=split(line, )
        PutSNPsIntoMatrix
        return
    return
"""
#Define functions
"""
This function reads a VCF file, remove SNPs of abnormal density (which are exported into a VCF file "OldFileName_removed.vcf"), 
and write other SNPs into a new SNP Matrix file "OldFileName_reserved.vcf". 

vcf_file_path: the file name of the input VCF file.
reserved_vcf_file_path: the file name with full path of the output VCF file containing all reserved SNPs.
removed_vcf_file_path: the file name with full path of the output VCF file containing all removed SNPs.
end_length: the length of the end regions in a contig, in which all SNPs will be removed. Default is 500.
window_size: the length of the window in which the number of SNPs should be no more than max_num_snp. Default is 1000.
max_num_snp: the maximum number of SNPs in a window. Default is 2.
"""
def FilterDenseSNPsFromVCF(vcf_file_path, reserved_vcf_file_path, removed_vcf_file_path, end_length=500, window_size=1000, max_num_snp=2):
    
    #reserved_vcf_file_path: the file name of the new VCF file containing all reserved SNPs. File name: OldVCFFileName_Reserved.vcf
    #reserved_vcf_file_path=vcf_file_path[:-4]+"_reserved.vcf"
    
    #removed_vcf_file_path: the file name of the VCF file containing all removed SNPs. File name: OldVCFFileName_Removed.vcf
    #removed_vcf_file_path=vcf_file_path[:-4]+"_removed.vcf"
    
    try:
        vcf_reader=vcf.Reader(open(vcf_file_path, 'r'))
    except:
        print "Cannot open the input vcf file.\n"
        return
    
    try:
        vcf_writer_reserved=vcf.Writer(open(reserved_vcf_file_path, 'w'), vcf_reader)
    except:
        print "Cannot create the file for reserved SNPs.\n"
        return
    
    try:
        vcf_writer_removed=vcf.Writer(open(removed_vcf_file_path, 'w'), vcf_reader)
    except:
        print "Cannot create the file for removed SNPs.\n"
        #close vcf_writer_reserved and remove the file reserved_vcf_file_path
        vcf_writer_reserved.close()
        os.remove(reserved_vcf_file_path)
        return
    
    #SNP list, saved as (Contig_Name, [(SNP_Position, SNP_Record),]), where SNP_Record is a line in VCF.
    snp_dict=dict()
            
    for vcf_data_line in vcf_reader:
        #get contig length from contig name.The CHROM should be a contig name in the format of Velvet/SPAdes output.
        try:
            ss=vcf_data_line.CHROM.split("_")
            contig_length=int(ss[ss.index("length")+1])
        except:
            #cannot find contig length. Should be a closed genome.
            contig_length=sys.maxsize
                
        if ((vcf_data_line.POS < end_length) or (vcf_data_line.POS > (contig_length - end_length))):
            #Position is close to the start of a contig
            #Add to removed_vcf_file_path with note "bad position"????
            vcf_writer_removed.write_record(vcf_data_line)
        else:
            #Only SNPs not in the end region are kept for further filtering
            key = vcf_data_line.CHROM
            if key not in snp_dict:
                record = [(vcf_data_line.POS, vcf_data_line)]
            else:
                record = snp_dict[key]
                temp_record=(vcf_data_line.POS, vcf_data_line)
                record.append(temp_record)
            snp_dict[key] = record
     
    #Process SNPs per contig        
    for contig, snp_list in snp_dict.items():
               
        #sort all SNPs in this contig by position
        sorted_list=sorted(snp_list, key=lambda SNPs: SNPs[0])
        
        #total number of SNPs
        num_of_snp=len(sorted_list)
        
        #remember removed SNPs by their positions
        removed_snp=list()
        
        #Process SNPs
        for idx, snp in enumerate(sorted_list):
            if ((idx + max_num_snp) < num_of_snp):
                if ((snp[0] + window_size) < sorted_list[idx+max_num_snp][0]):
                    if snp[0] not in removed_snp:
                        #not removed, so add to reserved list.
                        vcf_writer_reserved.write_record(snp[1])
                else:
                    for temp_idx in range(idx, idx+max_num_snp):
                        if sorted_list[temp_idx][0] not in removed_snp:
                            vcf_writer_removed.write_record(sorted_list[temp_idx][1])
                            removed_snp.append(sorted_list[temp_idx][0])
            else:
                #last max_num_snp SNPs
                if snp[0] not in removed_snp:
                    #not removed, so add to reserved list.
                    vcf_writer_reserved.write_record(snp[1])
     
    #Done with SNPs. Now close two VCF files and exit.    
    vcf_writer_reserved.close()
    vcf_writer_removed.close()       
        
"""Testing only                   
def main():
    #parse command line options
    filelist=list()
    path="/Users/yu.wang/work/snp_density/test"
    for filename in os.listdir(path):
        current_file=os.path.join(path, filename)
        filelist.append(current_file)
        
    for filename in filelist:
        FilterDenseSNPsFromVCF(filename)#, end_length, window_size, max_num_snp)
    
    sys.exit(0)
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Remove abnormally dense SNPs from the input VCF file, save the reserved SNPs into a new VCF file,
                                                    and save the removed SNPs into another VCF file.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('inputVCFFileName',     type=str, metavar='INPUT_VCF_FILE_NAME', help='File name of the input VCF file name with full path (or current path if full path is not provided) which must exist.')
    parser.add_argument('outputReservedVCFFileName',    type=str, metavar='OUTPUT_RESERVED_VCF_FILE_NAME', help='File name of the out VCF file name with full path (or current path if full path is not provided) containing all reserved SNPs.')
    parser.add_argument('outputRemovedVCFFileName',   type=str, metavar='OUTPUT_REMOVED_VCF_FILE_NAME', help='File name of the out VCF file name with full path (or current path if full path is not provided) containing all removed SNPs.')
    parser.add_argument('-l', '--end_length',  dest='endLength',    type=int, default=500, metavar='END_LENGTH', help='The length of the end regions in a contig, in which all SNPs will be removed. Default is 500.')
    parser.add_argument('-w', '--window_size',  dest='windowSize',  type=int, default=1000, metavar='WINDOW_SIZE', help='the length of the window in which the number of SNPs should be no more than max_num_snp. Default is 1000.')
    parser.add_argument('-m', '--max_snp',  dest='maxSNP',  type=int, default=2, metavar='MAX_NUM_SNPs', help='The maximum number of SNPs allowed in a window. Default is 2.')
    
    args=parser.parse_args()
    
    FilterDenseSNPsFromVCF(args.inputVCFFileName, args.outputReservedVCFFileName, args.outputRemovedVCFFileName, args.endLength, args.windowSize, args.maxSNP)
    sys.exit(0)