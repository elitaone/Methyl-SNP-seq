# python2.7
'''
Log:
Released on March 18, 2022
'''

'''
Requirements:
(1) Based on python2.7;
(2) bowtie2; command 'bowtie2' should be executable in $PATH, otherwise need to speficy the path.
(3) The following python scripts should be in the working directory (e.g. python FastqHandle_pairend.py), otherwise need to specify the path:
    FastqHandle_pairend.py, CompareBase.py, BaseCalibration.py, DeconvolutionCalibration_v2.py

Purpose:
Run DeconvolutionPipeline.py to script bash scripts to deconvolute the Methyl-SNP-seq Reads (Deconvolution Step) with Calibration:
(1) convert the bisulfited converted Non-methyl C base in Read1 back to C if R1=T and R2=C
    Set the base quality of R1 bases that are mismatching the R2 bases to a low base quality score based on the comparison to reference genome. 
(2) save cytosine methylation status in a report
(3) generate a fastq file Deconvolution_R1.fq and a methylation report

DeconvolutionPipeline.py to create a bash script DeconvolutionPipelinePrefix, including the following steps.
(1) Generate BaseCalibration table and probability table based on Deconvolution_R1 and R2
    Also perform removal of PCR duplicates and selection of unimapping here. Could be skipped to speed up.
(2) Deconvolute Read1 and Read2 including TtoC conversion, methylation calling and Base calibration.
    Output fastq file name_Deconvolution_Calibration_R1.fq is used for bowtie2 mapping. 
    Output methylation report name.Deconvolution.5mC is used for methylation calling.

Usage:
Run the following command to create bash scripts, then run the generated bash script:
$python DeconvolutionWithCalibration --Read1 TestSeq_hairpin_R1_val_1.fq --Read2 TestSeq_hairpin_R2_val_2.fq --name TestSeq -–reference hg38.fa

--Read1, Read2:
    illumina adapter and hairpin adapter removed Read1 and Read2 fastq files
--percent: Default 0.05
    percent of downsample of fastq for Base Calibration analysis
--reference: 
    reference genome.fa used to map the reads for calibration
--name: name of output files
    output files are saved in the current working directory:

    name.BaseCalibration.table
    name.BaseCalibration.probability
    name.Deconvolution.5mC: methylation report, cytosine position in read is 0-indexed.
    name_DeconvolutedRead.fq: Deconvolution for Read1 with BaseCalibration
--vcf: Optional
    A vcf file containing known SNP positions.
    If provided, the positions in this vcf file are ignored from base calibration.
--smp: Optional
    number of threads used for bowtie2 mapping, default 1
--dir: Optional
    directory to save the output files
    if not provided, the output files are saved at current working directory.
--path_to_python, --path_to_bowtie2:
    use this option to specify a path to the python2.7 and bowtie2 executable, e.g.
    /usr/bin/python and /usr/bin/bowtie2
    Else it is assumed that python2.7 and bowtie2 is executable in the PATH.
'''
import time
import os, sys
import argparse
from subprocess import check_call

script_dir = os.path.dirname( __file__ )
sys.path.append(os.path.join(script_dir, 'src'))

def Arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--name', help='name suffix used for output', dest='name')
    parser.add_argument('--percent', help='fastq downsample percent used for BaseCalibration analysis', dest='percent', default=0.05)
    parser.add_argument('--reference', help='reference fa file for mapping', dest='reference') 
    parser.add_argument('--vcf', help='vcf file showing snp', dest='vcf', required=False, default = None)
    parser.add_argument('--dir', help='dir to save the output file', dest='dir', required=False, default = None)
    parser.add_argument('--smp', help='number of thread used for bowtie2 mapping', dest='smp', type=int, default = 1)

    parser.add_argument('--path_to_python', help='specify a path to python2.7', dest='pathPython', default='python')
    parser.add_argument('--path_to_bowtie2', help='specify a path to bowtie2', dest='pathBowtie2', default='bowtie2')

    args = parser.parse_args()
    return args

def Deconvolution(Read1, Read2, percent, name, smp, dir, prefix, reference, vcf, pathPython, pathBowtie2):
    '''
    Generate a script DeconvolutionPipelinePrefix for several deconvolution steps.

    Read1, Read2: Read1 and Read2 fastq files in which the illumina adapter and hairpin adapter are removed
    percent: percent used to downsample of trimmred paired end fastq for BaseCalibration
    name: e.g. AccRep1
    dir: dir to save the output files
    '''
    
    with open('DeconvolutionPipeline{}'.format(prefix), 'w') as output:
        print "Generate script {} for Deconvolution with calibrarion.".format('DeconvolutionPipeline{}'.format(prefix))
        print>>output, '#!/bin/bash'
        
        CWD = os.getcwd()
        command = 'toucht DeconvolutionDir{}'.format(prefix)
        print>>output, command
        command = 'cd DeconvolutionDir{}'.format(prefix)
        print>>output, command

        # choose a portion of input reads pairs for calibration
        command = '{} FastqHandle_pairend.py --no-gzip downsample \
--input {} {} --output {}.downsample.R1.fq {}.downsample.R2.fq --percent {}'.format(pathPython, Read1, Read2, name, name, percent)
        print>>output, command
        # name.downsample.R1.fq, name.downsample.R2.fq

        # Convertion of Read1 to generate Deconvolution_R1
        command = '{} DeconvolutionConversion.py --Read1 {}.downsample.R1.fq --Read2 {}.downsample.R2.fq --name {}'.format(pathPython, name, name, name) 
        print>>output, command
        # name_Deconvolution_R1.fq, note this is named for downsampled Deconvolution_R1

        # bowtie2 mapping for downsampled Deconvolution_R1 and R2
        command = '{} {} {}_reference'.format(pathBowtie2.replace('bowtie2', 'bowtie2-build') ,reference, name) # bowtie2-build reference name
        print>>output, command
        command = '{} -p {} -x  {}_reference -U {}_Deconvolution_R1.fq -S {}_deconvoluted_R1.sam'.format(pathBowtie2, smp, name, name, name) 
        print>>output, command
        # name_deconvoluted_R1.sam
        command = '{} -p {} -x {}_reference -U {}.downsample.R2.fq -S {}_R2.sam'.format(pathBowtie2, smp, name, name, name) 
        print>>output, command
        # name_R2.sam
        
        # Compare the Deconvolution_R1, Read2 base with the reference for base calibration
        command = '{} CompareBase.py --bam {}_deconvoluted_R1.sam \
--fastq {}_Deconvolution_R1.fq --output {}_deconvoluted_R1.compareBase.txt --vcf {}'.format(pathPython, name, name, name, vcf)
        # name_deconvoluted_R1.compareBase.txt
        print>>output, command

        command = '{} CompareBase.py --bam {}_R2.uni.nodup.sam \
--fastq {}.downsample.R2.fq --output {}_R2.compareBase.txt --vcf {}'.format(pathPython, name, name, name, vcf)
        print>>output, command
        # name_R2.compareBase.txt

        command = '{} BaseCalibration.py \
--Read1 {}_deconvoluted_R1.compareBase.txt --Read2 {}_R2.compareBase.txt --name {}'.format(pathPython, name, name, name) # output: name.BaseCalibration.table, name.BaseCalibration.probability
        print>>output, command

        # Deconvolution including Base Calibration
        command = '{} DeconvolutionCalibration_v2.py --Read1 {} --Read2 {} --name {} \
--probability {}.BaseCalibration.probability'.format(pathPython, Read1, Read2, name, name) # name_Deconvolution_Calibration_R1.fq, name.Deconvolution.5mC
        print>>output, command
        command = 'mv {}_Deconvolution_Calibration_R1.fq {}_DeconvolutedRead.fq'.format(name, name)
        print>>output, command
        

        # move the output files to dir
        print>>output, 'echo The output files are saved at: {}'.format(dir)
        for item in ['{}_DeconvolutedRead.fq'.format(name), '{}.Deconvolution.5mC'.format(name), '{}.BaseCalibration.table'.format(name), '{}.BaseCalibration.probability'.format(name)]:
            command = 'mv {} {}'.format(item, dir)
            print>>output, command
        command = 'cd {}'.format(CWD)
        print>>output, command
        command = 'rm -r DeconvolutionDir{}'.format(prefix)
        print>>output, command

    return 'DeconvolutionPipeline{}'.format(prefix)
        

    
##-------mainbody
if __name__ == '__main__':

    # check the required python files:
    path = os.path.join(script_dir, 'src')
    pythonls = ['FastqHandle_pairend.py', 'DeconvolutionConversion.py', 'CompareBase.py', 'BaseCalibration.py', 'DeconvolutionCalibration_v2.py']

    for item in pythonls:
        if not os.path.exists(os.path.join(path, item)):
            print "{} can not be found.".format(item)
            quit()

    args = Arg()

    Read1 = os.path.abspath(args.Read1)
    Read2 = os.path.abspath(args.Read2)

    CWD = os.getcwd()
    if args.dir:
        dir = args.dir
    else:
        dir = CWD

    localtime = time.asctime(time.localtime())
    prefix = ''.join(localtime.split()[-2].split(':')) # '151542'

    script = Deconvolution(Read1, Read2, args.percent, args.name, args.smp, dir, prefix, args.reference, args.vcf, args.pathPython, args.pathBowtie2)
    check_call('bash {}'.format(script), shell=True)
