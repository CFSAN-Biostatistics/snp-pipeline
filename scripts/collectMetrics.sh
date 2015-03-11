#!/bin/bash
#$ -pe mpi 8
#$ -N extract_data_for_pt_report.sh
#$ -o /scratch/tmp/pt_data.out
#$ -e /scratch/tmp/pt_data.err
#$ -cwd

# Author: Joseph D. Baugher
# Copyright (c) 2015 Joseph D. Baugher (<joseph.baugher@fda.hhs.gov>). 
# All rights reserved.


#-------------------------
# Default Variables
#-------------------------

lab='Unknown'
output_file=''


#------
# Usage
#------

display_usage() { 
  echo -e "\nThis script is designed to extract desired data for a specific PT sequencing run." 
  echo -e "Please provide the path to the run analysis directory which includes fastq.gz files and"
  echo -e "output files of the SNP pipeline."
  echo -e "\nUsage: $0 [options] <run_path>\n"
  echo -e "Command line switches are optional, although the laboratory name (-l) should be provided."
  echo -e "  -c  Path to configuration file."
  echo -e "      Default: /isilon/.gnome/gnome3/genomeTrakr_files/PT_2014/analyses/scripts/process_PT_runs.config."
  echo -e "  -d  Debug" 
  echo -e "  -l  The name of the laboratory which performed the sequencing."
  echo -e "  -o  The name and path of the output file."
  echo -e "      Defaults to STDOUT."
  echo -e "  -h  Displays this help message\n"
  exit 0;
} 

#--------
# Options
#--------

while getopts ":c:dl:o:h" opt; do
  case $opt in
    c) #set option "c"
      config_path=$OPTARG
      ;; 
    d) #set option "d"
      debug=1
      ;;
    l) #set option "l"
      lab=$OPTARG
      ;; 
    o) #set option "o"
      output_file=$OPTARG
      ;; 
    h) #help
      display_usage
      ;;  
   \?) #unrecognized option - show help
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      display_usage
      ;;
  esac
done

#----------
# Arguments
#----------

shift $((OPTIND-1))

run_path=${1%/}

if [ -z $run_path ]; then display_usage; fi
if [ -z $config_path ]; then 
  config_path='./scripts/process_PT_runs.config'
fi
if [ ! -e $config_path ]; then
  echo "A configuration file must be provided.\n"
  display_usage
fi

# Attach config file
. $config_path

run_name=${run_path##*/}
sample_name=${run_name%%_*}


#-------------------------
# Get machine and flowcell from fastq header
#-------------------------

# Check if fastq or fastq.gz exist and are populated
fastq_gz=(${run_path}/*.fastq.gz)
if [ -s ${fastq_gz[0]} ]; then   # Check if files are *.fastq.gz 
  fq_head=$(gzip -cd ${fastq_gz[0]} | head -1)
else
  echo "Populated fastq.gz files were not located." 1>&2 ; exit 1;
fi  

# @M01402:25:000000000-A6P1T:1:1101:11489:1986 1:N:0:2
machine=${fq_head#@}
machine=${machine%%:*}
fields=(${fq_head//:/ })
flowcell=${fields[2]}

#-------------------------
# Sum file sizes (MB) of paired fastq files 
#------------------------- 
for file in ${fastq_gz[@]}; do   
  size=$(du -h ${file} | cut -f1)
  size=${size%%M*}
  file_size=$((file_size+size))
done

#-------------------------
# Predict which organism was sequenced
#------------------------- 
if [ ! -e ${run_path}/taxonomy_results.txt ]; then  
  module load kraken/0.10.5-beta/kraken
  kraken --preload --db $TAX_DB --threads 8 --paired ${fastq_gz[0]} ${fastq_gz[1]} \
    > ${run_path}/taxonomy_results.txt
fi
if [ ! -e ${run_path}/probable_taxa.txt ]; then 
  probable_taxa.pl ${run_path}/taxonomy_results.txt \
    > ${run_path}/probable_taxa.txt
fi
probable_taxa=$(tail -1 ${run_path}/probable_taxa.txt | cut -f1)

#-------------------------
# Get metadata from annotation files 
#-------------------------
biosample=$(awk -F'\t' '{if ($1 == s) print $3 }' s="${sample_name}" \
  ${META_PATH}/${SAMPLE_META})
taxonomy=$(awk -F'\t' '{if ($1 == s) print $2 }' s="${sample_name}" \
  ${META_PATH}/${SAMPLE_META})

#-------------------------
# Calculate number of reads and %mapped from sam file
#-------------------------
module load samtools
nreads=$(samtools view -c ${run_path}/${SAM})
mapped=$(samtools view -c -F 4 ${run_path}/${SAM})
perc_mapped=$(bc <<< "scale=4; ($mapped/$nreads)*100")
perc_mapped=$(printf '%.02f' $perc_mapped)

#-------------------------
# Count number of SNPs from vcf file
#-------------------------
snps=$(grep -v '^#' ${run_path}/${VCF} | wc -l)

#-------------------------
# Calculate mean depth from pileup file
#-------------------------    
depth=$(awk '{sum+=$4}END{print sum/NR}' ${run_path}/${PILEUP})
depth=$(printf '%.0f' $depth)

#-------------------------
# Get sequencing date from bionumerics dump (tab delim)
#-------------------------   
# Correct alternate newline characters 
seq_date=$(tr "\r" "\n" < ${BIONUMERICS} | tr -d '\015' | grep ${run_name} \
  | cut -f57)
# If multiple entries in bionumerics for same runid, choose the first one
seq_date=(${seq_date//\n/ })
seq_date=${seq_date[0]}
if [ -z "$seq_date" ]; then seq_date="NA"; fi

#-------------------------
# Get sample type from run.meta file, assumed to be 'Culture' unless otherwise specified.
#-------------------------   
sample_type=$(grep ${run_name} ${META_PATH}/${RUN_META} | cut -f2)
if [ -z "$sample_type" ]; then sample_type="Culture"; fi

#-------------------------
# Print Results
#-------------------------   
if [ -z $output_file ]; then 
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$lab" "$biosample" "$taxonomy" "$run_name" "$sample_type" "$file_size" \
  "$machine" "$flowcell" "$seq_date" "$nreads" "$perc_mapped" "$depth" "$snps" \
  "$probable_taxa"
else
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$lab" "$biosample" "$taxonomy" "$run_name" "$sample_type" "$file_size" \
  "$machine" "$flowcell" "$seq_date" "$nreads" "$perc_mapped" "$depth" "$snps" \
  "$probable_taxa" > $output_file
fi
  

debug() {
  echo "lab $lab"
  echo "biosample $biosample"
  echo "expected taxonomy $taxonomy"
  echo "run_name $run_name"
  echo "sample_type $sample_type"
  echo "file_size $file_size"
  echo "machine $machine"
  echo "flowcell $flowcell"
  echo "seq_date $seq_date"
  echo "nreads $nreads"
  echo "perc_mapped $perc_mapped"
  echo "depth $depth"
  echo "snps $snps"
  echo "probable_taxa $probable_taxa"
}

if [[ $debug ]]; then debug; fi    