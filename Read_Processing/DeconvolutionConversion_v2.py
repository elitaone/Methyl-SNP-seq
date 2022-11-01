# -*- coding: utf-8 -*-
# python2.7
'''
Made on Jan 26, 2022
'''
'''
Requirement: Based on python2.7

For methyl-SNP-Seq, run Read Deconvolution without Calibration.

Purpose:
For Methyl-SNP-seq paired end reads, compare positions in Read1 and Read2 to:
(1) convert the bisulfited converted Non-methyl C base in Read1 back to C if R1=T and R2=C
    Set the base quality of R1 bases that are mismatching the R2 bases to 0 (!). 
(2) save cytosine methylation status in a report
This step generates a fastq file name_DeconvolutedRead.fq and a methylation report name.Deconvolution.5mC
    
Logic:
    Criteria to report a read in output:
        match R1 and R2 from the first bp to the end, if len(R1)!=len(R2), the end/read length is the min(len(R1), len(R2))
        total R1, R2 mismatch < 10% of read length (100 base sequencing having 10 mismatches)
        here mismatch is defined as: R1!=R2 but except for R1==T and R2==C
    
    Convert a read by comparing R1 with R2:
    If R1!=N and R2!=N:
        if R1==R2:
            DeconvolutionBase=R1, DeconvolutionBaseQuality=Quality_R1
            if R1==R2==C, this position is marked as a methylated C in report
        else R1!=R2:
            if R1==T and R2==C:
                DeconvolutionBase=C, DeconvolutionBaseQuality=Quality_R1, this position is marked as an unmethylated C in report
            elif R1==C and R2!=C:
                DeconvolutionBase=R1=C, DeconvolutionBaseQuality=0
                this position is marked as a methylated C
            else:
                other cases R1!=R2 (e.g. R1==G and R2==A)
                DeconvolutionBase=R1, DeconvolutionBaseQuality=0
        else:
            if R1==N and R2!=N:
                DeconvolutionBase=R2, DeconvolutionBaseQuality=Quality_R2
                if R2==C, position is marked as an unmethylated C in report 
            else R1!=N and R2==N:
                DeconvolutionBase=R1, DeconvolutionBaseQuality=Quality_R1
                if R1==C, this position is marked as a methylated C

Usage:
$python DeconvolutionConversion_v2.py --Read1 TestSeq_R1_val_1.fq --Read2 TestSeq_R2_val_2.fq --name TestSeq

--Read1 and --Read2:
    fastq files for Read1 fastq and Read2 fastq in which the illumina adapter and hairpin adapter are removed.
    Only take unconpressed the fastq files.

--name: name of output files, e.g. TestSeq
    output files are saved in the current working directory.
    TestSeq_DeconvolutedRead.fq: Deconvoluted Read1 with Non-methy C convertd back and sequencing quality calibrated for mismatch bases
    TestSeq.Deconvolution.5mC: C or 5mC report

Note:
(1) Current version does not use multiple processing.
'''

try:
    import argparse
    import os
    import re
except:
    print "Import module error."
    quit()

__version__ = '2022.08.20'

class Deconvolute:
    '''
    Convert Read1 T to C if R1=T and R2=C
    '''
    def __init__(self, Read1, Read2, name):
        self.Read1 = os.path.abspath(Read1) # hairpin adapter removed Read1
        self.Read2 = os.path.abspath(Read2) # hairpin adapter removed Read2
        self.name = name

    def DeconvoluteEachLine(self, Read1, Read2, Read1_qual, Read2_qual):
        '''
        Use to deconvolute by comparing Read1 with Read2, and Recalibrate the sequencing quality for bases with R1, R2 mismatch
        len(Deconvolution Read0 = min(len(Read1), len(Read2))
        See above Logic for Deconvolution Rules.
        '''
        ls_R1=[] # list saving deconvolution_R1
        ls_R1_qual=[] # list saving quality score corresponding to deconvoluted base
        ls_C=[] # list saving C or Methy-C information
        count_mismatch=0

        for i in range(0, min(len(Read1), len(Read2))):
            if Read1[i] != 'N' and Read2[i] != 'N': 
                if Read1[i] == Read2[i]:
                    ls_R1.append(Read1[i])
                    ls_R1_qual.append(Read1_qual[i])
                    if Read1[i] == 'C':
                        ls_C.append('M{}'.format(i))
                else:
                    if Read1[i] == 'T' and Read2[i] == 'C':
                        ls_R1.append('C')
                        ls_R1_qual.append(Read1_qual[i])
                        ls_C.append('C{}'.format(i))
                    elif Read1[i] == 'C' and Read2[i] != 'C':
                        ls_R1.append(Read1[i])
                        ls_R1_qual.append('0') # DeconvolutionBaseQuality=0
                        ls_C.append('M{}'.format(i))
                        count_mismatch += 1
                    else: # other R1!=R2
                        ls_R1.append(Read1[i])
                        ls_R1_qual.append('0') # DeconvolutionBaseQuality=0
                        count_mismatch += 1
            else:
                if Read1[i] == 'N': # Modify N using R2 base for R1
                    ls_R1.append(Read2[i])
                    ls_R1_qual.append(Read2_qual[i])
                    if Read2[i] == 'C': # R1=N, R2=C, labeled as Non-methy C
                        ls_C.append('C{}'.format(i))
                else: # Do not modify N in R2
                    ls_R1.append(Read1[i])
                    ls_R1_qual.append(Read1_qual[i])
                    if Read1[i] == 'C': # R1=C, R2=N, labeled as methy C
                        ls_C.append('M{}'.format(i))

        return ''.join(ls_R1), ''.join(ls_R1_qual), ''.join(ls_C), count_mismatch, min(len(Read1), len(Read2))
        # deconvoluted_R1, deconvoluted_R1_quality, methylation information, count_mismatch, count_matching_length

    def DeconvoluteFastq(self):
        '''
        If count_mismatch/count_matching_length >= 10% (10bp mismatch in 100bp pair), exlude this unmatched pair
        Do not apply filter on Read having certain number of methy C.
        '''
        output_R1 = open('{}_Deconvolution_R1.fq'.format(self.name), 'w')
        output_methylation = open('{}.Deconvolution.5mC'.format(self.name), 'w')
        
        total_input = 0 # total number of read pairs in input
        total_base_mismatch = 0 # total bases of R1=!R2, excluding R1=T, R2=C and N
        total_base_matching_length = 0 # total bases for R1 used to compare to R2, length of match pair, excluding unmatched pair
        total_more_than_5mC = 0 # total number of reads having more than 5mC
        total_unmatched = 0 # total number of unmatched paired reads
        total_match = 0 # total number of match paired reads
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
                total_input += 1
                Decon_R1, Decon_R1_qual, methylation, count_mismatch, count_matching_length = \
                self.DeconvoluteEachLine(R1_line2, R2_line2, R1_line4, R2_line4)
                if float(count_mismatch)/count_matching_length < 0.1: # mismatch_percent_cutoff=0.1 could be changed
                    print>>output_R1, R1_line1
                    print>>output_R1, Decon_R1
                    print>>output_R1, R1_line3
                    print>>output_R1, Decon_R1_qual
                    print>>output_methylation, '\t'.join([R1_line1.split()[0].replace('@', ''), methylation])
                    if len(re.findall('M', methylation)) > 5:
                        total_more_than_5mC += 1
                    total_match += 1
                    total_base_matching_length += count_matching_length
                    total_base_mismatch += count_mismatch
                else:
                    total_unmatched += 1

                R1_line1 = f1.readline().strip()
                R1_line2 = f1.readline().strip()
                R1_line3 = f1.readline().strip()
                R1_line4 = f1.readline().strip()

                R2_line1 = f2.readline().strip()
                R2_line2 = f2.readline().strip()
                R2_line3 = f2.readline().strip()
                R2_line4 = f2.readline().strip()

            output_R1.close()
            output_methylation.close()
        
        print 'Perform Deconvolution without calibrarion using verion: {}.'.format(__version__)
        print "Read Pairs in input: {}".format(total_input)
        print "Matched Deconvoluted Read Pairs: {}".format(total_match)
        print "Unmatched Read Pairs: {}".format(total_unmatched)
        print "Read having more than 5 5mC: {}".format(total_more_than_5mC)
        print "Deconvoluted bases: {}".format(total_base_matching_length)
        print "Bases of R1!=R2: {}".format(total_base_mismatch)

        return 1

    
##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--name', help='name suffix used for output', dest='name')
    args = parser.parse_args()

    temp = Deconvolute(args.Read1, args.Read2, args.name)
    temp.DeconvoluteFastq()