#!/usr/bin/env python

from __future__ import print_function
import sys
import unittest
import argparse
import filecmp
import os
import subprocess
from testfixtures import TempDirectory
from pkg_resources import resource_filename
from snppipeline import cfsan_snp_pipeline

data_directory = resource_filename(cfsan_snp_pipeline.__name__, 'data')

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
        with open(in_path, "r") as in_file:
            for line in in_file:
                found = False
                for s in not_match_list:
                    if s in line:
                        found = True
                        break
                if not found:
                    out_file.write(line)


class SnpPipelineTest(unittest.TestCase):
    '''Unit test for snppipeline.'''

    def setUpOnce(self):
        """Create the list of sample directories"""
        print("Python version = " + sys.version)
        samples_directory = os.path.join(self.directory_run_result, 'samples')
        subdir_list = [os.path.join(samples_directory, subdir) for subdir in os.listdir(samples_directory)]

        with open(self.file_of_directories, "w") as path_file_object:
            for subdir in sorted(subdir_list):
                path_file_object.write("%s\n" % subdir)


    def compare_file(self, file_to_compare, not_match_str_list=None):
        """Compare a generated file to the expected results.
        """
        result_file = os.path.join(self.directory_run_result, file_to_compare)

        # Strip unwanted lines, like dates and versions, before the file compare
        correct_result_file = os.path.join(self.directory_correct, file_to_compare)
        self.assertTrue(os.path.isfile(correct_result_file), "The known correct result file does not exist: %s" % file_to_compare)
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


    def run_function_test(self, funct, args, file_to_compare, not_match_str_list=None):
        """Test one function and compare the generated file to the expected results.
        """
        args.verbose = 0

        # remember the previous file timestamp, if any
        result_file = os.path.join(self.directory_run_result, file_to_compare)
        pre_exist = os.path.isfile(result_file)
        if pre_exist:
            old_timestamp = os.path.getmtime(result_file)

        # Run the function under test and verify the output exists
        funct(args)
        self.assertTrue(os.path.isfile(result_file), "Result file does not exist: %s" % file_to_compare)

        # Verify result file has a newer timestamp if it previously existed
        if pre_exist:
            new_timestamp = os.path.getmtime(result_file)
            self.assertTrue(new_timestamp > old_timestamp, "Old file was not overwritten with newer timestamp: %s" % file_to_compare)

        # Strip unwanted lines, like dates and versions, and compare file
        self.compare_file(file_to_compare, not_match_str_list)


class SnpPipelineLambdaVirusTest(SnpPipelineTest):
    '''Unit test for snppipeline run with Lambda Virus data.'''

    all_setup_done = False

    def setUpOnce(self):
        """Create directories and data files for subsequent tests.
        """
        print("\nPreparing data files for tests.  This will take a minute...")
        temp_dir = TempDirectory()
        lambda_dir = os.path.join(temp_dir.path, "testLambdaVirus")
        ret = subprocess.call(["cfsan_snp_pipeline", "data", "lambdaVirusInputs", lambda_dir])
        SnpPipelineTest.directory_correct = os.path.join(lambda_dir, "lambdaVirusExpectedResults")
        ret = subprocess.call(["cfsan_snp_pipeline", "data", "lambdaVirusExpectedResults", SnpPipelineTest.directory_correct])
        samples_dir = os.path.join(lambda_dir, "samples")
        reference_file = os.path.join(lambda_dir, "reference", "lambda_virus.fasta")

        devNull = open(os.devnull, 'w')
        ret = subprocess.call("run_snp_pipeline.sh -o %s -s %s %s" % (lambda_dir, samples_dir, reference_file), shell=True, stdout=devNull)
        devNull.close()
        SnpPipelineTest.directory_run_result = lambda_dir
        SnpPipelineTest.file_of_directories = os.path.join(lambda_dir, 'sampleDirectories.txt')

        # Let the super class do some one-time setup
        SnpPipelineTest.setUpOnce(self)
        SnpPipelineLambdaVirusTest.all_setup_done = True



    @staticmethod
    def tearDownAll():
        """
        Delete all the temporary directories and files created during this
        testing session.
        """
        TempDirectory.cleanup_all()


    def __init__(self, *args):
        super(SnpPipelineLambdaVirusTest, self).__init__(*args)


    def setUp(self):
        """Setup"""
        if not SnpPipelineLambdaVirusTest.all_setup_done:
            self.setUpOnce()
            self.compare_file('metrics.tsv')


    def test_1_filter_regions(self):
        """Run filter_regions and verify var.flt_preserved.vcf and var.flt_removed.vcf contains expected contents for each sample.
        """
        command_line = "filter_regions -f " + \
                       os.path.join(self.__class__.directory_run_result, 'sampleDirectories.txt') + ' ' + \
                       os.path.join(self.__class__.directory_run_result, "reference", "lambda_virus.fasta")
        args = cfsan_snp_pipeline.parse_command_line(command_line)

        for dir in ['samples/sample1', 'samples/sample2','samples/sample3','samples/sample4']:
            self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, os.path.join(dir, 'var.flt_preserved.vcf'))
            self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, os.path.join(dir, 'var.flt_removed.vcf'))


    def test_2_merge_sites(self):
        """Run merge_sites and verify snplist.txt contains expected contents.
        """
        command_line = "merge_sites -f -o " + os.path.join(self.__class__.directory_run_result, 'snplist.txt') + ' ' + \
                       os.path.join(self.__class__.directory_run_result, 'sampleDirectories.txt') + ' ' + \
                       os.path.join(self.__class__.directory_run_result, 'filteredSampleDirectories.txt')
        args = cfsan_snp_pipeline.parse_command_line(command_line)
        self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, 'snplist.txt')


    def test_3_call_consensus(self):
        """Run call_consensus and verify consensus.fasta and consensus.vcf contain expected contents for each sample.
        """
        command_line = "call_consensus -f -l " + os.path.join(self.__class__.directory_run_result, 'snplist.txt') + \
                       " --vcfRefName lambda_virus.fasta" + \
                       " --vcfFailedSnpGt 1" + \
                       " dummyAllPileupFile"
        args = cfsan_snp_pipeline.parse_command_line(command_line)

        ignore_lines = ["##fileDate", "##source"]

        for dir in ['samples/sample1', 'samples/sample2','samples/sample3','samples/sample4']:
            args.allPileupFile = os.path.join(self.__class__.directory_run_result, dir, 'reads.all.pileup')
            args.consensusFile = os.path.join(self.__class__.directory_run_result, dir, 'consensus.fasta')
            args.vcfFileName = None
            self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, os.path.join(dir, 'consensus.fasta'))
            args.vcfFileName = 'consensus.vcf'
            self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, os.path.join(dir, 'consensus.vcf'), ignore_lines)


    def test_4_create_snp_matrix(self):
        """Run snp_matrix and verify snpma.fasta contains expected contents.
        """
        command_line = "snp_matrix -f -o " + os.path.join(self.__class__.directory_run_result, 'snpma.fasta') + ' ' + \
                        os.path.join(self.__class__.directory_run_result, 'sampleDirectories.txt')
        args = cfsan_snp_pipeline.parse_command_line(command_line)
        self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, 'snpma.fasta')


    def test_5_create_snp_reference_seq(self):
        """Run create_snp_reference_seq and verify referenceSNP.fasta contains expected contents.
        """
        command_line = "snp_reference -f -l " + os.path.join(self.__class__.directory_run_result, 'snplist.txt') + \
                       " -o " + os.path.join(self.__class__.directory_run_result, 'referenceSNP.fasta') + ' ' + \
                       os.path.join(self.__class__.directory_run_result, 'reference/lambda_virus.fasta')
        args = cfsan_snp_pipeline.parse_command_line(command_line)
        self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, 'referenceSNP.fasta')


    def test_6_calculate_snp_distances(self):
        """Run calculate_snp_distances and verify snp_distance_pairwise.tsv contains expected contents.
        """
        command_line = "distance -f " + os.path.join(self.__class__.directory_run_result, 'snpma.fasta') + ' ' + \
                       "--pairs " + os.path.join(self.__class__.directory_run_result, 'snp_distance_pairwise.tsv')
        args = cfsan_snp_pipeline.parse_command_line(command_line)
        self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, 'snp_distance_pairwise.tsv')


    def test_7_calculate_snp_distances(self):
        """Run calculate_snp_distances and verify snp_distance_matrix.tsv contains expected contents.
        """
        command_line = "distance -f " + os.path.join(self.__class__.directory_run_result, 'snpma.fasta') + ' ' + \
                       "--matrix " + os.path.join(self.__class__.directory_run_result, 'snp_distance_matrix.tsv')
        args = cfsan_snp_pipeline.parse_command_line(command_line)
        self.run_function_test(cfsan_snp_pipeline.run_command_from_args, args, 'snp_distance_matrix.tsv')


    def test_999(self):
        """Tear down"""
        self.tearDownAll()




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


