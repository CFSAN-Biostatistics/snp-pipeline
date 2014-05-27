========
Usage
========

The SNP Pipeline is run from the Unix command line.

Command Syntax
--------------
TODO: Document the main snp pipeline script with command line arguments


Step-by-Step Workflow
---------------------

Step 1 - Prep work::

    mkdir junk #make a directory to work in
    cd junk    #change to it

Step 2 - Get sequence data::

    #echo -e "SRR955145\nSRR955146\nERR178930\nERR178931\n" > prepInput
    echo -e "ERR178926\nERR178927\nERR178928\nERR178929\nERR178930\n" > prepInput
    cat prepInput | xargs -n 1 scripts/prepSequenceData.sh NC_011149

Step 3 - Set up reference sequence::

    scripts/prepReference.sh NC_011149

Step 4 - Set up sample sequence.
Note: This could be run in parallel using gnu parallel (workstation) or qsub (PBS on HPC)::

    cat prepInput | xargs -n 1 scripts/prepSamples.sh NC_011149
        
Step 5 - Run snp pipeline (samtools pileup in parallel and combine alignment and pileup to generate snp matrix)::

    ls -d -1 --color=never $PWD/samples/* > path.txt
    scripts/runsnppipeline.py -n 10 -d ~/snppipeline/test/junk/ -f path.txt -r reference/NC_011149.fasta -l snplist.txt -a snpma.fasta -i True
 
