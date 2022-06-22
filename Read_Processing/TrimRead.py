# python2.7
'''
Requirements:
(1) Based on python2.7;
(2) TrimGalore>=0.6.4, installed in PATH otherwise specify the path

Purpose:
Generate and run a bash script to trim the illumina adapter and hairpin adapter.
The resulting adapter removed Read1 and Read2 fastq files are used for Read Deconvolution.

Usage:
$python TrimRead.py --Read1 AccRep1.1.fastq.gz --Read2 AccRep1.2.fastq.gz --name AccRep1

--Read1, Read2:
    illumina sequencing read1, read2 fastq or fastq.gz
--name: Name for ourput fastq files

--path_to_cutadapt, --path_to_trimgalore:
    Use this option to specify a path to the Cutadapt executable, 
    e.g. /mnt/home/yan/exe/TrimGalore-0.6.4/trim_galore. 
    Else it is assumed that Cutadapt or trim_galore is in the PATH.

Output files:
    name_hairpin_R1_val_1.fq, name_hairpin_R2_val_2.fq
'''

import time
import os
import argparse
from subprocess import check_call

def Trim(Read1, Read2, name, prefix, pathCutadapt, pathTrimgalore):
    '''
    Generate a script TrimPrefix for several trimming steps
    trim galore to remove illumina adapter, hairpin adapter
    Read1, Read2: illumina sequencing Read1 and Read2
    Generate: 
    AccRep1_R1_val_1.fq
    AccRep1_R2_val_2.fq
    '''
    
    Read1 = os.path.abspath(Read1)
    Read2 = os.path.abspath(Read2)

    with open('Trim{}'.format(prefix), 'w') as output:
        # Remove illumina adapter
        print "Generate script {} for trimming the illumina and hairpin adapter.".format('Trim{}'.format(prefix))
        print>>output, '#!/bin/bash'

        if not pathCutadapt:
            command = '{} --paired {} {} --max_n 5 --trim-n --stringency 3 --illumina \
--basename {}_illumina --nextseq 20'.format(pathTrimgalore, Read1, Read2, name)
            print>>output, command
            # output: name_illumina_R1_val_1.fq.gz, name_illumina_R2_val_2.fq.gz

            # Remove hairpin adapter
            command = '{} --paired -a ATTATGATGATGATGATGAGTGTTAGGTTTGTTGTTGTTGTTGTGGT -a2 ACCACAACAACAACAACAAACCTAACACTCATCATCATCATCATAAT \
{}_illumina_R1_val_1.fq.gz {}_illumina_R2_val_2.fq.gz --dont_gzip \
--basename {}_hairpin'.format(pathTrimgalore, name, name, name)
            print>>output, command
            # output: name_hairpin_R1_val_1.fq, name_hairpin_R2_val_2.fq
        else:
            command = '{} --paired {} {} --max_n 5 --trim-n --stringency 3 --illumina \
--basename {}_illumina --nextseq 20 --path_to_cutadapt {}'.format(pathTrimgalore, Read1, Read2, name, pathCutadapt)
            print>>output, command

            command = '{} --paired -a ATTATGATGATGATGATGAGTGTTAGGTTTGTTGTTGTTGTTGTGGT -a2 ACCACAACAACAACAACAAACCTAACACTCATCATCATCATCATAAT \
{}_illumina_R1_val_1.fq.gz {}_illumina_R2_val_2.fq.gz --dont_gzip \
--basename {}_hairpin --path_to_cutadapt {}'.format(pathTrimgalore, name, name, name, pathCutadapt)
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
    
    parser.add_argument('--path_to_cutadapt', help='specify a path to the Cutadapt executable', dest='pathCutadapt', default='cutadapt')
    parser.add_argument('--path_to_trimgalore', help='specify a path to the trimgalore executable', dest='pathTrimgalore', default='trim_galore')

    

    args = parser.parse_args()

    localtime = time.asctime(time.localtime())
    prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

    script = Trim(args.Read1, args.Read2, args.name, prefix, args.pathCutadapt, args.pathTrimgalore)
    check_call('sh {}'.format(script), shell=True)