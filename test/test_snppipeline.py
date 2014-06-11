#!/usr/bin/env python2.7

import unittest

import filecmp
import os
import inspect

from snppipeline import snppipeline

# script directory
test_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# various directories of test files
test_lambda_virus_directory = os.path.join(test_directory, 'testLambdaVirus')
test_agona_mom_directory = os.path.join(test_directory, 'testAgonaMOM')
compare_lambda_virus_directory = os.path.join(test_directory, 'codeComparisonFiles', 'testLambdaVirus')
compare_agona_mom_virus_directory = os.path.join(test_directory, 'codeComparisonFiles', 'testAgonaMOM')

class Test(unittest.TestCase):
    '''Unit test for snppipeline.'''

    def test_snppipeline_lambda_virus(self):
        """Run snppipeline with synthetic virus example.
        """

        args_dict = {
            'maxThread':3,
            'mainPath': test_lambda_virus_directory,
            'Reference':'lambda_virus.fa',
            'pathFileName':'path.txt',
            'snplistFileName':'snplist.txt',
            'snpmaFileName':'snpma.fasta',
            'bamFileName':'reads.bam',
            'pileupFileName':'reads.pileup',
            'verbose':False,
            'includeReference':True,
            'useOldPileups':False,
            'combinedDepthAcrossSamples':10,
            'alleleFrequencyForFirstALTAllele':1.0,
            'arFlagValue':1.0
        }

        snppipeline.run_snp_pipeline(args_dict)

        #Compare the files in the two directories whose names are given.
        #Returns three lists of file names: match, mismatch, errors.
        directory_correct = compare_lambda_virus_directory
        directory_run_result = test_lambda_virus_directory
        files_to_compare = ['snplist.txt','snpma.fasta','sample1/reads.pileup',
                            'sample2/reads.pileup','sample3/reads.pileup',
                            'sample4/reads.pileup','referenceSNP.fasta']
        match, mismatch, errors =  filecmp.cmpfiles(directory_correct,
                                                    directory_run_result,
                                                    files_to_compare,shallow=False)
        print('Match, Mismatch, Errors: '+str(len(match))+', '+str(len(mismatch))+', '+str(len(errors)))
        print('  Match: ',match)
        print('  Mismatch: ',mismatch)
        print('  Errors: ',errors)
        self.assertEqual(True,len(match)    == len(files_to_compare) and
                              len(mismatch) == 0 and
                              len(errors)   == 0)

        #Remove files generated #TODO make this optional?
        for file_name in files_to_compare:
            os.remove(os.path.join(directory_run_result,file_name))


    def test_snppipeline_agona(self):
        """Run snppipeline with agona 5 samples example.
        """

        args_dict = {
            'maxThread':8,
            'mainPath': test_agona_mom_directory,
            'Reference':'NC_011149.fasta',
            'pathFileName':'path.txt',
            'snplistFileName':'snplist.txt',
            'snpmaFileName':'snpma.fasta',
            'bamFileName':'reads.bam',
            'pileupFileName':'reads.pileup',
            'verbose':False,
            'includeReference':True,
            'useOldPileups':False,
            'combinedDepthAcrossSamples':10,
            'alleleFrequencyForFirstALTAllele':1.0,
            'arFlagValue':1.0
        }

        snppipeline.run_snp_pipeline(args_dict)

        #Compare the files in the two directories whose names are given.
        #Returns three lists of file names: match, mismatch, errors.
        directory_correct = compare_agona_mom_virus_directory
        directory_run_result = test_agona_mom_directory
        files_to_compare = ['snplist.txt',
                            'snpma.fasta',
                            'samples/CFSAN_genomes/CFSAN000448/reads.pileup',
                            'samples/CFSAN_genomes/CFSAN000449/reads.pileup',
                            'samples/CFSAN_genomes/CFSAN000450/reads.pileup',
                            'samples/SRA_data/ERR178930/reads.pileup',
                            'samples/SRA_data/ERR178931/reads.pileup',
                            'referenceSNP.fasta']
        match, mismatch, errors =  filecmp.cmpfiles(directory_correct,
                                                    directory_run_result,
                                                    files_to_compare,shallow=False)
        print('  Match: ',match)
        print('  Mismatch: ',mismatch)
        print('  Errors: ',errors)
        print('Match, Mismatch, Errors: '+str(len(match))+', '+str(len(mismatch))+', '+str(len(errors)))

        self.assertEqual(True,len(match)    == len(files_to_compare) and
                              len(mismatch) == 0 and
                              len(errors)   == 0)

        #TODO make this optional?
        for file_name in files_to_compare:
            os.remove(os.path.join(directory_run_result,file_name))

if __name__ == "__main__":
    unittest.main()


    #compare archived files to expected files

