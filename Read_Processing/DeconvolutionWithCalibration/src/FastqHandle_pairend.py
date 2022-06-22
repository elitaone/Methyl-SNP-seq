#!/mnt/home/ettwiller/yan/exe/miniconda2/envs/my-r-env/bin/python
# By Bo Yan
'''
Log:

Created on Aug 16, 2020 based on FastqHandle.py
Modifications on Nov 14, 2020:
    Add option --no-gzip to control the format of output fastq.
Modifications on Feb 9, 2022:
    Fix the bugs only the tempoary input files are removed in EndStep.
'''

"""
Based on Python2.7

This file contains several functions for handling the paired fastq files.

# Select read pairs with both Read1 and Read2 length >= cutoff to output file
$python FastqHandle_pairend.py select --input Read1 Read2 --output Read1.output Read2.output --cutoff int

# Downsample fastq of Read1 and Read2 to a given proportion together, so output files are still paired: output Read1 and Read2 having the same readID.
$python FastqHandle_pairend.py downsample --input Read1 Read2 --output Read1.output Read2.output --percent 0.35 float
$python FastqHandle_pairend.py --no-gzip downsample --input Read1 Read2 --output Read1.output Read2.output --percent 0.35 float
Roungly randomly select about 35% reads from input fastq, does not output a exact number.

--input: paired Read1 and Read2, fastq or fastq.gz
--output: 
    paired fastq or fastq.gz corresponding to Read1 and Read2
    if the input is gz, the output is gziped. The name should also contain .gz.
    if the input is fastq, the output is fastq too.
--no-gzip:
    Add --no-gzip to create the output files without compression.
    without this option, the output file will be compressed.
"""
try:
    import os
    import argparse
    from subprocess import check_call
    import time
    import random
except:
    print "module error."
    quit()

class Fastq:
    def __init__(self, input, output):
        localtime = time.asctime(time.localtime())
        self.suffix = ''.join(localtime.split()[-2].split(':'))
        
        self.input = [os.path.abspath(item) for item in input] # a list [Read1, Read2]
        self.input_fq, self.gz = self.FormatInputs() # self.input_fq is a tuple (Read1, Read2)

        self.output_temp = ["{}.read1.output".format(self.suffix), "{}.read2.output".format(self.suffix)] # temporary output fastq, e.g. 192306.read1.output, 192306.read2.output
        self.output = [os.path.abspath(item) for item in output] # a list [Read1.output, Read2.output] for output, gzip the fastq if input if compressed

    def CheckFormat(self, fileformat):
        '''
        check the format of fileformat
        '''
        # automatically detect gzip or not
        command = 'gzip -l {} > /dev/null'.format(fileformat)
        try:
            check_call(command, shell=True)
            print "Input {} is a fastq.gz file.".format(fileformat)
            return True # file is compressed
        except:
            return False

    def FormatInputs(self):
        '''
        uncompress the input Read1 and Read2 if they are compressed
        If Read1 and Read2 are in different compressed format, quit
        return ((Read1.temp.fq, Read2.temp.fq), True/False for gz)
        '''
        if all(self.CheckFormat(item) for item in self.input): # Both Read1 and Read2 are gz
            command = 'gunzip -c {} > {}.R1.fastq'.format(self.input[0], self.suffix) # uncompress the gz -> suffix.fastq
            check_call(command, shell=True)
            command = 'gunzip -c {} > {}.R2.fastq'.format(self.input[1], self.suffix) # uncompress the gz -> suffix.fastq
            check_call(command, shell=True)
            return (('{}.R1.fastq'.format(self.suffix), '{}.R2.fastq'.format(self.suffix)), True)
        elif all(not self.CheckFormat(item) for item in self.input): # Both Read1 and Read2 are not gz
            return ((self.input[0], self.input[1]), False)
        else:
            print "Read1 and Read2 are in different compressed format. Quit"
            quit()

    def Select(self, cutoff):
        '''
        only save the read pairs with length>=cutoff in both Read1 and Read2 in the output.
        '''
        print "Select the reads with length above {}.".format(cutoff)

        total_output = 0
        total_input = 0
        output1 = open(self.output_temp[0], 'w')
        output2 = open(self.output_temp[1], 'w')
        with open(self.input_fq[0]) as f1:
            with open(self.input_fq[1]) as f2:
                R1_line1 = f1.readline()
                R1_line2 = f1.readline()
                R1_line3 = f1.readline()
                R1_line4 = f1.readline()

                R2_line1 = f2.readline()
                R2_line2 = f2.readline()
                R2_line3 = f2.readline()
                R2_line4 = f2.readline()
                while R1_line1:
                    total_input += 1
                    if len(R1_line2.strip())>=cutoff and len(R2_line2.strip())>=cutoff:
                        print>>output1, R1_line1.strip()
                        print>>output1, R1_line2.strip()
                        print>>output1, R1_line3.strip()
                        print>>output1, R1_line4.strip()

                        print>>output2, R2_line1.strip()
                        print>>output2, R2_line2.strip()
                        print>>output2, R2_line3.strip()
                        print>>output2, R2_line4.strip()
                        
                        total_output += 1
                    R1_line1 = f1.readline()
                    R1_line2 = f1.readline()
                    R1_line3 = f1.readline()
                    R1_line4 = f1.readline()

                    R2_line1 = f2.readline()
                    R2_line2 = f2.readline()
                    R2_line3 = f2.readline()
                    R2_line4 = f2.readline()
        output1.close()
        output2.close()
        
        print "There are {} number of reads in the input Read1 or Read2.".format(total_input)        
        print "There are {} number of reads in the output Read1 or Read2.".format(total_output)


    def Downsample(self, percent):
        '''
        Downsample Read1 and Read2 together, to make sure the output files having the same and same order of readID
        '''
        percent_int = round(percent * 100) # return the nearest integer
        
        total_input = 0
        total_output = 0
        output1 = open(self.output_temp[0], 'w')
        output2 = open(self.output_temp[1], 'w')
        with open(self.input_fq[0]) as f1:
            with open(self.input_fq[1]) as f2:
                R1_line1 = f1.readline()
                R1_line2 = f1.readline()
                R1_line3 = f1.readline()
                R1_line4 = f1.readline()

                R2_line1 = f2.readline()
                R2_line2 = f2.readline()
                R2_line3 = f2.readline()
                R2_line4 = f2.readline()
                while R1_line1:
                    total_input +=1
                    if random.randint(1, 100) <= percent_int:
                        print>>output1, R1_line1.strip()
                        print>>output1, R1_line2.strip()
                        print>>output1, R1_line3.strip()
                        print>>output1, R1_line4.strip()

                        print>>output2, R2_line1.strip()
                        print>>output2, R2_line2.strip()
                        print>>output2, R2_line3.strip()
                        print>>output2, R2_line4.strip()

                        total_output +=1
                    R1_line1 = f1.readline()
                    R1_line2 = f1.readline()
                    R1_line3 = f1.readline()
                    R1_line4 = f1.readline()

                    R2_line1 = f2.readline()
                    R2_line2 = f2.readline()
                    R2_line3 = f2.readline()
                    R2_line4 = f2.readline()
        output1.close()
        output2.close()
        print "There are {} number of reads in the input Read1 or Read2.".format(total_input)
        print "There are {} number of reads in the output Read1 or Read2.".format(total_output)


    def EndStep(self, nogzip):
        if not nogzip: # compress the output in the absence of option --no-gzip
            print "The output file is compressed."
            command = 'gzip -c {} > {}'.format(self.output_temp[0], self.output[0])
            check_call(command, shell=True)
            command = 'gzip -c {} > {}'.format(self.output_temp[1], self.output[1])
            check_call(command, shell=True)
            os.remove(self.output_temp[0])
            os.remove(self.output_temp[1])
        else: # Do not compress the output in the presence of option --no-gzip
            print "The output file is not compressed."
            command = 'mv {} {}'.format(self.output_temp[0], self.output[0])
            check_call(command, shell=True)
            command = 'mv {} {}'.format(self.output_temp[1], self.output[1])
            check_call(command, shell=True)
        if self.gz: # input files are compressed
            os.remove(self.input_fq[0])
            os.remove(self.input_fq[1])
        print "The output files are: {} and {}".format(self.output[0], self.output[1])
    

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--no-gzip', help = 'Do not compress output file', dest = 'nogzip', action= 'store_true')

    subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')

    parser_a = subparsers.add_parser('select', help='select reads with length above cutoff')
    parser_a.add_argument('--input', nargs="+", help = 'input fastq.gz or fastq', dest='input_file')
    parser_a.add_argument('--output', nargs="+", help='output fastq.gz or fastq', dest='output_file')
    parser_a.add_argument('--cutoff', help='length cutoff', dest='cutoff', type=int)

    parser_b = subparsers.add_parser('downsample', help='downsample to a proportion')
    parser_b.add_argument('--input', nargs="+", help = 'input fastq.gz or fastq', dest='input_file')
    parser_b.add_argument('--output', nargs="+", help='output fastq.gz or fastq', dest='output_file')
    parser_b.add_argument('--percent', help='proportion for downsample', dest='percent', type=float)

    args = parser.parse_args()

    if len(args.input_file) < 2 or len(args.output_file) < 2:
        print "Input and Ouput need to have both Read1 and Read2."
        quit()

    if args.mode =='select':
        Fastqfile = Fastq(args.input_file, args.output_file)
        Fastqfile.Select(args.cutoff)
        Fastqfile.EndStep(args.nogzip)
    elif args.mode == 'downsample':
        Fastqfile = Fastq(args.input_file, args.output_file)
        Fastqfile.Downsample(args.percent)
        Fastqfile.EndStep(args.nogzip)


        

