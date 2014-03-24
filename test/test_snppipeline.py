#!/usr/bin/env python2.7

import unittest
#import snppipeline
import filecmp

import imp
snppipeline = imp.load_source('snppipeline', '/home/hugh.rand/projects/snppipeline/snppipeline/snppipeline.py')
utils       = imp.load_source('utils', '/home/hugh.rand/projects/snppipeline/snppipeline/utils.py')

class Test(unittest.TestCase):
    '''Unit tests for snppipeline.'''

    def test_snppipeline_lambda_virus(self):
        '''Run snppipeline with synthetic virus example.'''
        args_dict = {
            'maxThread':3,      
            'mainPath':'/home/hugh.rand/projects/snppipeline/test/testLambdaVirus/',        
            'Reference':'lambda_virus.fa',     
            'pathFileName':'path.txt',   
            'snplistFileName':'snplist.txt', 
            'snpmaFileName':'snpma.fa',
            'bamFileName':'reads.bam',
            'pileupFileName':'reads.pileup',
            'verbose':False,
            'includeReference':True,
            'useOldPileups':False
        } 
        snppipeline.run_snp_pipeline(args_dict)
        #TODO do diffs here
        #Compare the files in the two directories whose names are given.
        #Returns three lists of file names: match, mismatch, errors.
        directory_correct = '../test/codeComparisonFiles/testLambdaVirus' #TODO
        directory_run_result = '../test/testLambdaVirus' #TODO
        files_to_compare = ['snplist.txt','snpma.fasta','sample1/reads.pileup',
                            'sample2/reads.pileup','sample3/reads.pileup',
                            'sample4/reads.pileup','referenceSNP.fasta ']
        match, mismatch, errors =  filecmp.cmpfiles(directory_correct,
                                                    directory_run_result,
                                                    files_to_compare,shallow=False)
        print 'Match:', match
        print 'Mismatch:', mismatch
        print 'Errors:', errors
                                    
        self.assertEqual(True,match==7 and mismatch==0 and errors==0)
        pass

    def test_snppipeline_agona(self):
        '''Run snppipeline with synthetic virus example.'''
        args_dict = {
            'maxThread':3,      
            'mainPath':'~/projects/snppipeline/test/testLambdaVirus/',        
            'Reference':'reference.fasta',     
            'pathFileName':'path.txt',   
            'snplistFileName':'snplist.txt', 
            'snpmaFileName':'snpma.fa',
            'bamFileName':'reads.bam',
            'pileupFileName':'reads.pileup',
            'verbose':False,
            'includeReference':True,
            'useOldPileups':False
        } 
        #snppipeline.run_snp_pipeline(args_dict)
        pass

if __name__ == "__main__":
    unittest.main()
    

    #compare archived files to expected files

    #remove files from run    
    
#    def test_run_snp_pipeline(self)    
#
#class TestSequenceFunctions(unittest.TestCase):
#
#    def setUp(self):
#        self.seq = range(10)
#
#    def test_shuffle(self):
#        # make sure the shuffled sequence does not lose any elements
#        random.shuffle(self.seq)
#        self.seq.sort()
#        self.assertEqual(self.seq, range(10))
#
#        # should raise an exception for an immutable sequence
#        self.assertRaises(TypeError, random.shuffle, (1,2,3))
#
#    def test_choice(self):
#        element = random.choice(self.seq)
#        self.assertTrue(element in self.seq)
#
#    def test_sample(self):
#        with self.assertRaises(ValueError):
#            random.sample(self.seq, 20)
#        for element in random.sample(self.seq, 5):
#            self.assertTrue(element in self.seq)
#
#if __name__ == '__main__':
#    unittest.main()