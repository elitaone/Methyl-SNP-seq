# -*- coding: utf-8 -*-
# python2.7
'''
Made on Nov 7, 2020
Modified on Aug 20, 2022:
    change check_call to shell=False
    remove multiple threads since using downsampled reads
'''
'''
Based on python2.7

Purpose:
For Methyl-SNP-seq paired end reads, compare positions in Read1 and Read2.
This scrpit does not extract methylation information.
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
    output file is a fastq file: name_Deconvolution_R1.fq
'''

import argparse
import os

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
        output_R1 = open('{}_Deconvolution_R1.fq'.format(self.name), 'w')
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

    
##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--name', help='name suffix used for output', dest='name')
    args = parser.parse_args()

    Deconvolute(args.Read1, args.Read2, args.name)