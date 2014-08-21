#!/usr/bin/env python2.7

import unittest

import filecmp
import os
import inspect
from pkg_resources import resource_filename

from snppipeline import snppipeline

test_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
data_directory = resource_filename(snppipeline.__name__, 'data')

# various directories of test files
test_lambda_virus_directory = os.path.join(test_directory, 'testLambdaVirus')
test_agona_directory = os.path.join(test_directory, 'testAgona')
compare_lambda_virus_directory = os.path.join(data_directory, 'lambdaVirusExpectedResults')
compare_agona_directory = os.path.join(data_directory, 'agonaExpectedResults')
path_file_name = os.path.join(test_lambda_virus_directory, 'sampleDirectories.txt')

class SnpPipelineTest(unittest.TestCase):
    '''Unit test for snppipeline.'''

    def setUp(self):
        """Create the list of sample directories"""

        # TODO: Check to insure data files for the test exist.
        #       if not then tell user what script to run to create such a directory

        samples_directory = os.path.join(self.directory_run_result, 'samples')
        subdir_list = [os.path.join(samples_directory, subdir) for subdir in os.listdir(samples_directory)]

        with open(path_file_name, "w") as path_file_object:
            for subdir in sorted(subdir_list):
                path_file_object.write("%s\n" % subdir)


    def run_function_test(self, funct, args_dict, file_to_compare):
        """Test one function and compare the generated file to the expected results.
        """
        correct_result_file = os.path.join(self.directory_correct, file_to_compare)

        # remember the previous file timestamp, if any
        result_file = os.path.join(self.directory_run_result, file_to_compare)
        pre_exist = os.path.isfile(result_file)
        if pre_exist:
            old_timestamp = os.path.getmtime(result_file)

        funct(args_dict)
        self.assertTrue(os.path.isfile(result_file), "Result file does not exist: %s" % file_to_compare)
        match = filecmp.cmp(correct_result_file, result_file, shallow=False)
        self.assertTrue(match, "Incorrect file contents: %s" % file_to_compare)
        #print "\nMatch"
        #print correct_result_file
        #print result_file

        # Verify file has a newer timestamp if it previously existed
        if pre_exist:
            new_timestamp = os.path.getmtime(result_file)
            self.assertTrue(new_timestamp > old_timestamp, "Old file was not overwritten with newer timestamp: %s" % file_to_compare)


class SnpPipelineLambdaVirusTest(SnpPipelineTest):
    '''Unit test for snppipeline run with Lambda Virus data.'''

    def __init__(self, *args):
        self.directory_correct = compare_lambda_virus_directory
        self.directory_run_result = test_lambda_virus_directory
        super(SnpPipelineLambdaVirusTest, self).__init__(*args)

    def test_snppipeline_lambda_virus(self):
        """Run snppipeline with synthetic virus example.
        """

        args_dict = {
            'sampleDirsFile' : os.path.join(self.directory_run_result, 'sampleDirectories.txt'),
            'vcfFileName' : 'var.flt.vcf',
            'snpListFile' : os.path.join(self.directory_run_result, 'snplist.txt'),
            }
        self.run_function_test(snppipeline.create_snp_list, args_dict, 'snplist.txt')

        args_dict = {
            'snpListFile' : os.path.join(self.directory_run_result, 'snplist.txt'),
            }
        for dir in ['samples/sample1', 'samples/sample2','samples/sample3','samples/sample4']:
            args_dict['allPileupFile'] = os.path.join(self.directory_run_result, dir, 'reads.all.pileup')
            args_dict['snpPileupFile'] = os.path.join(self.directory_run_result, dir, 'reads.snp.pileup')
            self.run_function_test(snppipeline.create_snp_pileup, args_dict, os.path.join(dir, 'reads.snp.pileup'))
        
        args_dict = {
            'sampleDirsFile' : os.path.join(self.directory_run_result, 'sampleDirectories.txt'),
            'snpListFile' : os.path.join(self.directory_run_result, 'snplist.txt'),
            'pileupFileName' : 'reads.snp.pileup',
            'snpmaFile' : os.path.join(self.directory_run_result, 'snpma.fasta'),
            }
        self.run_function_test(snppipeline.create_snp_matrix, args_dict, 'snpma.fasta')

        args_dict = {
            'referenceFile' : os.path.join(self.directory_run_result, 'reference/lambda_virus.fasta'),
            'snpListFile' : os.path.join(self.directory_run_result, 'snplist.txt'),
            'snpRefFile' : os.path.join(self.directory_run_result, 'referenceSNP.fasta'),
            } 
        self.run_function_test(snppipeline.create_snp_reference_seq, args_dict, 'referenceSNP.fasta')




#TODO uncomment the agona test, and maybe make it optional depending on 
#       the availability of the data for running the test?
#
#class SnpPipelineAgonaTest(SnpPipelineTest):
#    '''Unit test for snppipeline run with Agona MOM data.'''
#
#    def __init__(self, *args):
#        self.directory_correct = compare_agona_directory
#        self.directory_run_result = test_agona_directory
#        super(SnpPipelineAgonaTest, self).__init__(*args)
#
#    def test_snppipeline_agona(self):
#        """Run snppipeline with agona 5 samples example.
#        """
#
#        args_dict = {
#            'maxThread':8,      
#            'mainPath': test_agona_directory,
#            'Reference':'NC_011149.fasta',     
#            'pathFileName':'sampleDirectoryNames.txt',   
#            'snplistFileName':'snplist.txt', 
#            'snpmaFileName':'snpma.fasta',
#            'bamFileName':'reads.bam',
#            'pileupFileName':'reads.pileup',
#            'verbose':False,
#            'includeReference':True,
#            'useOldPileups':False,
#        } 
#
#        #TODO Add test to insure data for test is in directory specified in args_dict['mainPath'],
#        #       if not then tell user what script to run to create such a directory
#
#        snppipeline.run_snp_pipeline(args_dict)
#        
#        #Compare the files in the two directories whose names are given.
#        #Returns three lists of file names: match, mismatch, errors.
#        self.directory_correct = compare_agona_virus_directory
#        self.directory_run_result = test_agona_directory
#        files_to_compare = ['snplist.txt',
#                            'snpma.fasta',
#                            'samples/CFSAN000448/reads.pileup',  
#                            'samples/CFSAN000449/reads.pileup', 
#                            'samples/CFSAN000450/reads.pileup',
#                            'samples/ERR178930/reads.pileup',
#                            'samples/ERR178931/reads.pileup',
#                            'referenceSNP.fasta']
#        match, mismatch, errors =  filecmp.cmpfiles(self.directory_correct,
#                                                    self.directory_run_result,
#                                                    files_to_compare,shallow=False)
#        print('  Match: ',match)
#        print('  Mismatch: ',mismatch)
#        print('  Errors: ',errors)
#        print('Match, Mismatch, Errors: '+str(len(match))+', '+str(len(mismatch))+', '+str(len(errors)))
#        
#        self.assertEqual(True,len(match)    == len(files_to_compare) and
#                              len(mismatch) == 0 and
#                              len(errors)   == 0)
#        
#        #TODO make this optional?
#        for file_name in files_to_compare:
#            os.remove(os.path.join(self.directory_run_result,file_name)) 

if __name__ == "__main__":
    unittest.main()


