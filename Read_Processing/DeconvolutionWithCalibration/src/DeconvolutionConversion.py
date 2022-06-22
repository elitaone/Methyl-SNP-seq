#!/mnt/home/ettwiller/yan/exe/miniconda2/envs/my-python2-env/python
'''
Made on Nov 7, 2020
'''
'''
Requirement:
Based on python2.7, need pathos

Purpose:
For Accurateseq paired end reads, compare positions in Read1 and Read2.
This is the same as Yi's hairpin_pair_directional_match_and_deconvolute_with_quality_5mC_v2020.pl,
but does not extract methylation information.
The generated Deconvolution_R1.fq could be used to calculate the BaseCalibration.probability for DeconvolutionCalibration.py.

Logic:

    Criteria for report read in Deconvolution read:
        match R1 and R2 from the first bp to the end, if len(R1)!=len(R2), the end/read length is the min(len(R1), len(R2))
        total R1, R2 mismatch < 10% of read length (100 base sequencing having 10 mismatches)
        here mismatch is defined as: R1!=R2 but except R1==T and R2==C
    
    Deconvolution Read by comparing R1 with R2:
    if R1==R2:
        DeconvolutionBase_R1=R1, DeconvolutionBaseQuality=Quality_R1
    else R1!=R2:
        if R1==T and R2==C:
            DeconvolutionBase_R1=C, DeconvolutionBaseQuality=Quality_R1
        else other R1!=R2 (e.g. R1==G and R2==A):
            DeconvolutionBase_R1=R1, DeconvolutionBaseQuality=Quality_R1

Usage:
$python /mnt/home/ettwiller/yan/Bo_script/DeconvolutionConversion.py \
--Read1 /home/yan/Bo/AccurateSeq/NA12878_AccurateSeq/deconvoluted/AccRep1_R1_val_1.fq \
--Read2 /home/yan/Bo/AccurateSeq/NA12878_AccurateSeq/deconvoluted/AccRep1_R2_val_2.fq \
--name AccRep1

--name: name of output files
    output file is a fastq file:
    name_Deconvolution_R1.fq

Note:
(1) Use multiple thread pathos.pools (nodes=ncpu) to speed up.
(2) tempory dir is named based on time, could be called by different runs.
(3) tempory dir is about 3 times large as the Read1 fastq, e.g. AccRep1_R1_trim.fq=250G, tempdir=700G. Watch out the space.
'''

try:
    import argparse
    import time
    import os
    from subprocess import check_call
    import re
    from pathos.pools import ProcessPool
except:
    print "Import module error."
    quit()

class Deconvolute:
    '''
    Convert Read1 T to C if R1=T and R2=C
    '''
    def __init__(self, Read1, Read2, name):
        self.Read1 = os.path.abspath(Read1) # hairpin adapter removed Read1
        self.Read2 = os.path.abspath(Read2) # hairpin adapter removed Read2
        self.name = name # name of the output Deconvlution files: Deconvolution_R1.name
        
        self.count = self.DeconvoluteFastq() # total number of deconvoluted reads in output

    def DeconvoluteEachLine(self, Read1, Read2, Read1_qual):
        '''
        Use to deconvolute by comparing Read1 with Read2
        See above Logic for Deconvolution Rules.
        '''
        ls_R1=[] # list saving deconvolution_R1
        ls_R1_qual=[] # list saving quality score corresponding to deconvoluted base
        count_mismatch=0

        for i in range(0, min(len(Read1), len(Read2))):
            if Read1[i] == Read2[i]:
                ls_R1.append(Read1[i])
            else:
                if Read1[i] == 'T' and Read2[i] == 'C':
                    ls_R1.append('C')
                else:
                    count_mismatch += 1
                    ls_R1.append(Read1[i])
            ls_R1_qual.append(Read1_qual[i])
            
        return ''.join(ls_R1), ''.join(ls_R1_qual), count_mismatch, min(len(Read1), len(Read2))
        # deconvoluted_R1, deconvoluted_R1_quality, count_mismatch, count_matching_length

    def DeconvoluteFastq(self):
        '''
        If count_mismatch/count_matching_length >= 10% (10bp mismatch in 100bp pair), exlude this unmatched pair
        Do not apply filter on Read having certain number of methy C.
        '''
        output_R1 = open('Deconvolution_R1.{}'.format(self.name), 'w')
        total_output = 0
        with open(self.Read1) as f1, open(self.Read2) as f2:
            R1_line1 = f1.readline().strip()
            R1_line2 = f1.readline().strip()
            R1_line3 = f1.readline().strip()
            R1_line4 = f1.readline().strip()
            
            R2_line1 = f2.readline().strip()
            R2_line2 = f2.readline().strip()
            R2_line3 = f2.readline().strip()
            R2_line4 = f2.readline().strip()
            while R1_line1:
                Decon_R1, Decon_R1_qual, count_mismatch, count_matching_length = \
                self.DeconvoluteEachLine(R1_line2, R2_line2, R1_line4)
                if float(count_mismatch)/count_matching_length < 0.1: # mismatch_percent_cutoff=0.1 could be changed
                    print>>output_R1, R1_line1
                    print>>output_R1, Decon_R1
                    print>>output_R1, R1_line3
                    print>>output_R1, Decon_R1_qual
                    total_output += 1
                R1_line1 = f1.readline().strip()
                R1_line2 = f1.readline().strip()
                R1_line3 = f1.readline().strip()
                R1_line4 = f1.readline().strip()

                R2_line1 = f2.readline().strip()
                R2_line2 = f2.readline().strip()
                R2_line3 = f2.readline().strip()
                R2_line4 = f2.readline().strip()

            output_R1.close()

        return total_output


def main(Read1, Read2, name):
    '''
    Split the large Read1 and Read2 file into file with 100 million reads
    Run class DeconvoluteCalibration for each using multiple threads
    '''
    localtime = time.asctime(time.localtime())
    prefix = ''.join(localtime.split()[-2].split(':')) # '151542'
    
    Read1 = os.path.abspath(Read1)
    Read2 = os.path.abspath(Read2)
    CWD = os.getcwd()
 
    # Using multiple thread
    os.mkdir('Deconvolution{}'.format(prefix)) # Deconvolution172545
    os.chdir('Deconvolution{}'.format(prefix))


    command = 'split -l 400000000 -d -a 2 {} R1_'.format(Read1) # R1_00, R1_01, R1_02 ...
    check_call(command, shell=True)
    command = 'split -l 400000000 -d -a 2 {} R2_'.format(Read2) # R2_00, R2_01, R2_02 ...
    check_call(command, shell=True)

    ls_input_R1 = []
    ls_input_R2 = []
    for item in os.listdir(os.getcwd()):
        if 'R1_' in item:
            ls_input_R1.append(item)
        elif 'R2_' in item:
            ls_input_R2.append(item)

    pool = ProcessPool() # nodes=ncpu
    ls_output_name = [item.replace('R1_', '') for item in sorted(ls_input_R1)] # [00, 01, 02 ...]
    ls_return = pool.map(Deconvolute, sorted(ls_input_R1), sorted(ls_input_R2), ls_output_name)
    # output files: Deconvolution_R1.00
    # ls_return: a list of call Deconvolute object, each has property count: total number of deconvoluted reads in output
   
    pool.close()
    pool.join()

    # combine all the output files 
    with open(os.path.join(CWD, '{}_Deconvolution_R1.fq'.format(name)), 'w') as output:
        for input in ['Deconvolution_R1.{}'.format(x) for x in ls_output_name]: # Deconvolution_R1 for split files
            with open(input) as f:
                for line in f:
                    print>>output, line.strip()

    total_output = sum([item.count for item in ls_return])

    print "Total number of deconvoluted Reads in output fastq: {}".format(total_output)

    os.chdir(CWD)
    command = 'rm -r {}'.format('Deconvolution{}'.format(prefix))
    check_call(command, shell=True)
    return 1

    
##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--name', help='name suffix used for output', dest='name')
    args = parser.parse_args()

    main(args.Read1, args.Read2, args.name)