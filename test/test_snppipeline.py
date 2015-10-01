#!/usr/bin/env python2.7

import unittest

import filecmp
import os
import subprocess
from testfixtures import TempDirectory
from pkg_resources import resource_filename
from snppipeline import snppipeline

data_directory = resource_filename(snppipeline.__name__, 'data')

# various directories of test files
compare_agona_directory = os.path.join(data_directory, 'agonaExpectedResults')


def grep_not_matching(in_path, out_path, not_match_list):
    """
    Copy all lines in a file except those containing any string in a given
    list of strings.

    Parameters
    ----------
    in_path : str
        Path to input file
    out_path : str
        Path to output file
    not_match_list : list of str
        List of strings to exclude from the copy
    """
    with open(out_path, "w") as out_file:
        for line in open(in_path):
            found = False
            for s in not_match_list:
                if s in line:
                    found = True
                    break
            if not found:            
                out_file.write(line)


class SnpPipelineTest(unittest.TestCase):
    '''Unit test for snppipeline.'''

    def setUp(self):
        """Create the list of sample directories"""
        samples_directory = os.path.join(self.directory_run_result, 'samples')
        subdir_list = [os.path.join(samples_directory, subdir) for subdir in os.listdir(samples_directory)]

        with open(self.file_of_directories, "w") as path_file_object:
            for subdir in sorted(subdir_list):
                path_file_object.write("%s\n" % subdir)


    def run_function_test(self, funct, args_dict, file_to_compare, not_match_str_list=None):
        """Test one function and compare the generated file to the expected results.
        """
        # remember the previous file timestamp, if any
        result_file = os.path.join(self.directory_run_result, file_to_compare)
        pre_exist = os.path.isfile(result_file)
        if pre_exist:
            old_timestamp = os.path.getmtime(result_file)

        # Run the function under test and verify the output exists
        funct(args_dict)
        self.assertTrue(os.path.isfile(result_file), "Result file does not exist: %s" % file_to_compare)

        # Verify result file has a newer timestamp if it previously existed
        if pre_exist:
            new_timestamp = os.path.getmtime(result_file)
            self.assertTrue(new_timestamp > old_timestamp, "Old file was not overwritten with newer timestamp: %s" % file_to_compare)

        # Strip unwanted lines, like dates and versions, before the file compare
        correct_result_file = os.path.join(self.directory_correct, file_to_compare)
        if not_match_str_list:
            correct_result_file2 = os.path.join(self.directory_correct, file_to_compare+"2")
            grep_not_matching(correct_result_file, correct_result_file2, not_match_str_list)
            correct_result_file = correct_result_file2

            result_file2 = os.path.join(self.directory_run_result, file_to_compare+"2")
            grep_not_matching(result_file, result_file2, not_match_str_list)
            result_file = result_file2

        # Compare files
        match = filecmp.cmp(correct_result_file, result_file, shallow=False)
        self.assertTrue(match, "Incorrect file contents: %s" % file_to_compare)
        #print "\nMatch"
        #print correct_result_file
        #print result_file


class SnpPipelineLambdaVirusTest(SnpPipelineTest):
    '''Unit test for snppipeline run with Lambda Virus data.'''

    @classmethod
    def setUpClass(cls):
        """Create directories and data files for subsequent tests.
        """

        print "Preparing data files for tests.  This will take a minute..."
        temp_dir = TempDirectory()
        lambda_dir = os.path.join(temp_dir.path, "testLambdaVirus")
        ret = subprocess.call(["copy_snppipeline_data.py", "lambdaVirusInputs", lambda_dir])
        cls.directory_correct = os.path.join(lambda_dir, "lambdaVirusExpectedResults")
        ret = subprocess.call(["copy_snppipeline_data.py", "lambdaVirusExpectedResults", cls.directory_correct])
        samples_dir = os.path.join(lambda_dir, "samples")
        reference_file = os.path.join(lambda_dir, "reference", "lambda_virus.fasta")

        devNull = open(os.devnull, 'w')
        ret = subprocess.call("run_snp_pipeline.sh -o %s -s %s %s" % (lambda_dir, samples_dir, reference_file), shell=True, stdout=devNull)
        devNull.close()
        cls.directory_run_result = lambda_dir
        cls.file_of_directories = os.path.join(lambda_dir, 'sampleDirectories.txt')


    @classmethod
    def tearDownClass(self):
        """
        Delete all the temporary directories and files created during this 
        testing session.
        """
        TempDirectory.cleanup_all()


    def __init__(self, *args):
        super(SnpPipelineLambdaVirusTest, self).__init__(*args)


    def test_1_create_snp_list(self):
        """Run create_snp_list and verify snplist.txt contains expected contents.
        """
        args_dict = {
            'sampleDirsFile' : os.path.join(self.__class__.directory_run_result, 'sampleDirectories.txt'),
            'vcfFileName' : 'var.flt.vcf',
            'snpListFile' : os.path.join(self.__class__.directory_run_result, 'snplist.txt'),
            'forceFlag' : True,
            }
        self.run_function_test(snppipeline.create_snp_list, args_dict, 'snplist.txt')


    def test_2a_create_snp_pileup(self):
        """Run create_snp_pileup and verify reads.snp.pileup contains expected contents for each sample.
        """
        args_dict = {
            'snpListFile' : os.path.join(self.__class__.directory_run_result, 'snplist.txt'),
            'forceFlag' : True,
            }
        for dir in ['samples/sample1', 'samples/sample2','samples/sample3','samples/sample4']:
            args_dict['allPileupFile'] = os.path.join(self.__class__.directory_run_result, dir, 'reads.all.pileup')
            args_dict['snpPileupFile'] = os.path.join(self.__class__.directory_run_result, dir, 'reads.snp.pileup')
            self.run_function_test(snppipeline.create_snp_pileup, args_dict, os.path.join(dir, 'reads.snp.pileup'))


    def test_2b_call_consensus(self):
        """Run call_consensus and verify consensus.fasta and consensus.vcf contain expected contents for each sample.
        """
        args_dict = {
            'snpListFile' : os.path.join(self.__class__.directory_run_result, 'snplist.txt'),
            'forceFlag' : True,
            }

        ignore_lines = ["##fileDate", "##source"]

        for dir in ['samples/sample1', 'samples/sample2','samples/sample3','samples/sample4']:
            args_dict['allPileupFile'] = os.path.join(self.__class__.directory_run_result, dir, 'reads.all.pileup')
            args_dict['consensusFile'] = os.path.join(self.__class__.directory_run_result, dir, 'consensus.fasta')
            args_dict['minBaseQual'] = 0
            args_dict['minConsFreq'] = 0.6
            args_dict['minConsStrdDpth'] = 0
            args_dict['minConsStrdBias'] = 0
            args_dict['forceFlag'] = True
            args_dict['vcfFileName'] = None
            args_dict['vcfAllPos'] = False
            self.run_function_test(snppipeline.call_consensus, args_dict, os.path.join(dir, 'consensus.fasta'))
            args_dict['vcfFileName'] = 'consensus.vcf'
            args_dict['vcfRefName'] = 'lambda_virus.fasta'
            self.run_function_test(snppipeline.call_consensus, args_dict, os.path.join(dir, 'consensus.vcf'), ignore_lines)


    def test_3_create_snp_matrix(self):
        """Run create_snp_matrix and verify snpma.fasta contains expected contents.
        """
        args_dict = {
            'sampleDirsFile' : os.path.join(self.__class__.directory_run_result, 'sampleDirectories.txt'),
            'consFileName' : 'consensus.fasta',
            'snpmaFile' : os.path.join(self.__class__.directory_run_result, 'snpma.fasta'),
            'forceFlag' : True,
            }
        self.run_function_test(snppipeline.create_snp_matrix, args_dict, 'snpma.fasta')


    def test_4_create_snp_reference_seq(self):
        """Run create_snp_reference_seq and verify referenceSNP.fasta contains expected contents for each sample.
        """
        args_dict = {
            'referenceFile' : os.path.join(self.__class__.directory_run_result, 'reference/lambda_virus.fasta'),
            'snpListFile' : os.path.join(self.__class__.directory_run_result, 'snplist.txt'),
            'snpRefFile' : os.path.join(self.__class__.directory_run_result, 'referenceSNP.fasta'),
            'forceFlag' : True,
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


