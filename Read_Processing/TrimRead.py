# -*- coding: utf-8 -*-
# python2.7
'''
Created on June 7, 2022
Modifications on Nov 1, 2022:
    Add step to get input file path and confirm the existence.
'''
'''
Requirements:
(1) Based on python2.7;
(2) TrimGalore>=0.6.4, installed in PATH otherwise specify the path

Purpose:
Generate and run a bash script to trim the illumina adapter and hairpin adapter.
The resulting adapter removed Read1 and Read2 fastq files are used for Read Deconvolution.

Usage:
$python TrimRead.py --Read1 TestSeq.1.fastq.gz --Read2 TestSeq.2.fastq.gz --name TestSeq

--Read1, Read2:
    illumina sequencing read1, read2 fastq or fastq.gz

--name: name for ourput fastq files, e.g. TestSeq
Output files: TestSeq_hairpin_R1_val_1.fq, TestSeq_hairpin_R2_val_2.fq

--cores: Default 1
    Number of cores to be used for trim_galore trimming.
    It is used as –cores/-j option for trim_galore.

--path_to_cutadapt, --path_to_trimgalore:
    Use this option to specify a path to the trim_galore or Cutadapt executable, 
    e.g. /mnt/home/yan/exe/TrimGalore-0.6.4/trim_galore. 
    Else it is assumed that Cutadapt or trim_galore is in the PATH.

Note:
If this pipeline does not work well (likely due to the naming or path error), you can also run trim_galore by yourself in two steps as shown below:
Step1 to remove illumina adapter:
$trim_galore --paired Read1.fastq.gz Read2.fastq.gz --max_n 5 --trim-n --stringency 3 --illumina --nextseq 20 --cores 2
The trimming parameter can be adjusted based on your sequencing.

Step2 to remove hairpin adapter:
$trim_galore --paired -a ATTATGATGATGATGATGAGTGTTAGGTTTGTTGTTGTTGTTGTGGT -a2 ACCACAACAACAACAACAAACCTAACACTCATCATCATCATCATAAT \
R1_val_1.fq.gz R2_val_2.fq.gz --dont_gzip --cores 2 –-basename NameYouWant
Here R1_val_1.fq.gz R2_val_2.fq.gz are the output files of Step1 trimming.
'''

import time
import os
import argparse
from subprocess import check_call

__version__ = '2022.11.01'

def GetFilePath(input_file):
    '''
    return the abspath of the input_file
    if use relative path and input_file is not in cwd, the abspath will be wrong, return False
    '''
    filepath = os.path.abspath(input_file)
    
    if os.path.exists(filepath):
        return filepath
    else:
        print('Can not find file: {}'.format(filepath))
        return None

def Trim(Read1, Read2, name, prefix, pathCutadapt, pathTrimgalore, cores):
    '''
    Generate a script TrimPrefix for several trimming steps
    trim galore to remove illumina adapter, hairpin adapter
    Read1, Read2: illumina sequencing Read1 and Read2
    Generate: 
    AccRep1_R1_val_1.fq
    AccRep1_R2_val_2.fq
    '''
    
    Read1 = GetFilePath(Read1)
    Read2 = GetFilePath(Read2)

    with open('Trim{}'.format(prefix), 'w') as output:
        # Remove illumina adapter
        print "Generate script {} for trimming the illumina and hairpin adapter.".format('Trim{}'.format(prefix))
        print>>output, '#!/bin/bash'

        if not pathCutadapt:
            # Remove the illumina adapters
            command = '{} --paired {} {} --max_n 5 --trim-n --stringency 3 --illumina \
--basename {}_illumina --nextseq 20 --cores {}'.format(pathTrimgalore, Read1, Read2, name, cores)
            print>>output, command
            # output: name_illumina_R1_val_1.fq.gz, name_illumina_R2_val_2.fq.gz

            # Remove the hairpin adapter
            command = '{} --paired -a ATTATGATGATGATGATGAGTGTTAGGTTTGTTGTTGTTGTTGTGGT -a2 ACCACAACAACAACAACAAACCTAACACTCATCATCATCATCATAAT \
{}_illumina_R1_val_1.fq.gz {}_illumina_R2_val_2.fq.gz --dont_gzip \
--basename {}_hairpin --cores {}'.format(pathTrimgalore, name, name, name, cores)
            print>>output, command
            # output: name_hairpin_R1_val_1.fq, name_hairpin_R2_val_2.fq
        else:
            command = '{} --paired {} {} --max_n 5 --trim-n --stringency 3 --illumina \
--basename {}_illumina --nextseq 20 --path_to_cutadapt {} --cores {}'.format(pathTrimgalore, Read1, Read2, name, pathCutadapt, cores)
            print>>output, command

            command = '{} --paired -a ATTATGATGATGATGATGAGTGTTAGGTTTGTTGTTGTTGTTGTGGT -a2 ACCACAACAACAACAACAAACCTAACACTCATCATCATCATCATAAT \
{}_illumina_R1_val_1.fq.gz {}_illumina_R2_val_2.fq.gz --dont_gzip \
--basename {}_hairpin --path_to_cutadapt {} --cores {}'.format(pathTrimgalore, name, name, name, pathCutadapt, cores)
            print>>output, command

        print "The output files {} and {} are saved under {}".format('{}_hairpin_R1_val_1.fq'.format(name), '{}_hairpin_R2_val_2.fq'.format(name), os.getcwd())

        command = 'rm {}_illumina_R1_val_1.fq.gz'.format(name)
        print>>output, command
        command = 'rm {}_illumina_R2_val_2.fq.gz'.format(name)
        print>>output, command

        command = 'rm {}'.format('Trim{}'.format(prefix))
        print>>output, command

    return 'Trim{}'.format(prefix)

##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--name', help='name suffix used for output fastq file', dest='name')
    parser.add_argument('--cores', help='Number of cores to be used for trimming', dest='cores', default=1, type=int)
    
    parser.add_argument('--path_to_cutadapt', help='specify a path to the Cutadapt executable', dest='pathCutadapt', default='cutadapt')
    parser.add_argument('--path_to_trimgalore', help='specify a path to the trimgalore executable', dest='pathTrimgalore', default='trim_galore')

    args = parser.parse_args()

    localtime = time.asctime(time.localtime())
    prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

    script = Trim(args.Read1, args.Read2, args.name, prefix, args.pathCutadapt, args.pathTrimgalore, args.cores)
    check_call(['sh', script], shell=False)
