# -*- coding: utf-8 -*-
# python2.7
'''
# Logs:

Created on June 16, 2021
'''
'''
Requirement:
Provide the path of samtools in SortBamFile function.

Explanation:
This script is used to correct some of the unmatched read pairs that Read1 and Read2 bases are not matched properly due to insertion of deletion.
The output name_unmatchcorrect.R1.fastq, name_unmatchcorrect.R1.fastq files can be deconvoluted using DeconvolutionWithCalibration or DeconvolutionConversion_v2.py.

Logic: 
(1) extract unmatched read pairs that can not be deconvoluted by DeconvolutionWithCalibration.
    unmatched R1 and R2 can be extract using CompareTwoFastq.py, e.g.
    $python CompareTwoFastq.py different --bait TarRep_400bp_Deconvolution_Calibration_R1.fq \
    --target TarRep_400bp_R1_trim.fq --output TarRep_400bp_R1_unmatch.fq
    unmatched R1 and R2 are mapped using bismark and bowtie2, respectively.
(2) unmatch R1, R2 read pairs that are qualified for correction if:
    R1, R2 are mapped to the same chr and in the same direction
    there are pairs that R1 and R2 are Reverse complementary, which map to the same loci in opposite direction
    The distance between R1 mapping start and R2 mapping start is <=5
(3) Pair the SEQ_Index to REF_Position based on the mapping of R1 and R2, use REF_Position to realign the bases in read pair
    If a mismatch of R1/R2 is an insertion relative to REF, remove this base and the corresponding qual from read
    If a mismatch of R1/R2 is a deletion in one of the read, use the base on the other read and the corresponding qual
    If a mismatch of R1/R2 is a deletion in both reads, pass

Usage exampple:
$python DeconvolutionUnmatchCorrect.py --Read1 TestSeqUnmatch_R1.bam --Read2 TestSeqUnmatch_R2.bam --name TestSeq

--Read1:
    unmatched R1 reads aligned using bismark mapping, the resulting sam or bam file
--Read2:    
    unmatched R2 reads aligned using bowtie2 mapping, the resulting sam or bam file
    The selection of -F 4 reads and sorting by name will be performed in the function SortBamFile in this script, 
    so do no need to pre selection or sorting for Read1 or Read2 file.
--name:
    Output fastq files are: name_unmatchcorrect.R1.fastq, name_unmatchcorrect.R2.fastq
--dir: Optional
    Dir for saving output files, not required. If not provided, output files are save in the current working dir.
--path_to_samtools: 
    Use this option to specify a path to the samtools executable,
    e.g. /usr/bin/samtools. Else by default it is assumed that samtools and bedtools is in the PATH.
    
Note:
(1) Examples of read pairs before and after correction.
e.g.
A00336:A00336:HF2VMDRXY:1:1101:10438:13839
Read1: TATAAAGAGTTTAGTTTTATATGTTTTTTTTTTTTT-ATTTTTTGTATTTATTATGAGTATAAGGTGATTTTCATTTTTAATTATTATTTGATGTTTAT
QUAL1: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFF,FFFFFFFFFFFFFFFFFF,FFFFFFFFFF:FFFFF
Read2: TACAAAGAGTTCAGCTTCATATGCTCCTTCCCCTTC ATTTCTTGCATTTATTATGAACATAAGGTGATTTTCATCTCCAATGATTACCTGATGCTTAT
QUAL2: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

After correction:
R1: TATAAAGAGTTTAGTTTTATATGTTTTTTTTTTTTTATTTTTTGTATTTATTATGAGTATAAGGTGATTTTCATTTTTAATTATTATTTGATGTTTAT
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFF,FFFFFFFFFFFFFFFFFF,FFFFFFFFFF:FFFFF
R2: TACAAAGAGTTCAGCTTCATATGCTCCTTCCCCTTCATTTCTTGCATTTATTATGAACATAAGGTGATTTTCATCTCCAATGATTACCTGATGCTTAT
'''

try:
    from subprocess import check_call
    import os
    import re
    import time
    import numpy as np
    import argparse
except:
    print "module error"
    quit()

__version__ ='2022.08.20'

def Flag16(Flag):
    '''
    return True if Flag contains Flag 16, Flase if not
    Flag 0 return False
    '''
    if len(str(np.binary_repr(int(Flag))))<5: # Flag < 16
        return False
    else:
        return str(np.binary_repr(int(Flag)))[-5] == '1'

def RC_seq(seq):
    '''
    reverse complementary conversion
    '''
    dic={'A':'T', 'T':'A', 'G':'C','C':'G', 'U':'T', 'N':'N'}
    seq = seq.upper()
    ls = list(seq)
    ls_new = []
    for letter in ls[::-1]:
        if dic.has_key(letter):
            ls_new.append(dic[letter])
        else:
            ls_new=[]
            print "Input seq contains non DNA."
            exit()
    return ''.join(ls_new)

# --------
def SortBamFile(Read1_bam, Read2_bam, suffix, pathSamtools):
    '''
    Extract the mappable reads (-F 4 -F 256) and sort the reads by name
    Sorting by name guarantee the reads in R1_forcorrect.sam and R2_forcorrect.sam have the same order

    By default, bismark output does not contain -f 4 reads but bowtie2 output does
    Generate and Return:
    R1_sort.suffix.sam, R2_sort.suffix.sam
    '''
    print "Sort the input bam files."

    with open('samtools.{}.script'.format(suffix), 'w') as output:
        print>>output, '#!/bin/sh'
        print>>output, '{} view -h -F 4 -F 256 {} | {} sort -n -o {} -'.format(pathSamtools, Read1_bam, pathSamtools, 'R1_sort.{}.sam'.format(suffix))
        print>>output, '{} view -h -F 4 -F 256 {} | {} sort -n -o {} -'.format(pathSamtools, Read2_bam, pathSamtools, 'R2_sort.{}.sam'.format(suffix))
    
    try:    
        check_call(['sh', 'samtools.{}.script'.format(suffix)], shell=False)
        os.remove('samtools.{}.script'.format(suffix))
    except:
        print "Samtools error."
        quit()
    return 'R1_sort.{}.sam'.format(suffix), 'R2_sort.{}.sam'.format(suffix)

# --------
def ExtractReadsForCorrection(bismark_R1_sam, bowtie_R2_sam, suffix):
    '''
    bismark_R1_sam: R1 sorted by name (-F 4, -F 256) sam file
    bowtie_R2_sam: R2 sorted by name (-F 4, -F 256) sam file

    Use to extract the R1, R2 read pairs qualified for correction:
        R1, R2 are mapped to the same chr and in the same direction
        there are pairs that R1 and R2 are Reverse complementary, which map to the same loci in opposite direction
        The distance between R1 mapping start and R2 mapping start is <=5

    Generate and return: 
    R1_forcorrect.sam and R2_forcorrect.sam, do not have header, have the same order of reads
    '''
    print "Extract Read Pairs for Correction."

    # mappable Read2
    with open(bowtie_R2_sam) as f:
        key = [line.strip().split()[0] for line in f if not line.startswith('@')] # A00336:A00336:HF2VMDRXY:1:1101:10031:20055
        f.seek(0)
        value = [(line.strip().split()[2], line.strip().split()[3], Flag16(line.strip().split()[1])) for line in f if not line.startswith('@')]
        dic_R2 = dict(zip(key, value)) # {ID: (chr, start, Flag16(Flag))}

    # dic_need_correct, readID for reads: R1, R2 mapped to the same chr and distance<=5
    ls_ID_correct = [] # readID for reads: R1, R2 mapped to the same chr and distance<=5
    # mappable Read1
    with open(bismark_R1_sam) as f:
        for line in f:
            if not line.startswith('@'):
                key = line.strip().split()[0].split('_')[0] # A00336:A00336:HF2VMDRXY:1:1101:10004:12994_1:N:0:CTCGAACAGGAATTGC 
                if dic_R2.has_key(key):
                    # on the same chr, Flag16(FlagR1)==Flag16(FlagR2) return True if maping of R1 and R2 are in the same direction, False if they are in opposition direction
                    if line.strip().split()[2] == dic_R2[key][0] and Flag16(line.strip().split()[1])==dic_R2[key][2]:
                        distance = int(line.strip().split()[3]) - int(dic_R2[key][1])
                        if abs(distance)<=5:
                            ls_ID_correct.append(key)
    dic_ID_correct = dict.fromkeys(ls_ID_correct)

    # R1_forcorrect.sam and R2_forcorrect.sam do not have header, have the same order of reads
    count = 0
    output_R1 = open('R1_forcorrect.{}.sam'.format(suffix), 'w')
    with open(bismark_R1_sam) as f1:
        for line in f1:
            if not line.startswith('@') and dic_ID_correct.has_key(line.strip().split()[0].split('_')[0]):
                print>>output_R1, line.strip()
                count += 1
    output_R1.close()
    print "There are {} number of Read1-Read2 pairs are corrected for deconvolution conversion.".format(count)

    output_R2 = open('R2_forcorrect.{}.sam'.format(suffix), 'w')
    with open(bowtie_R2_sam) as f2:
        for line in f2:
            if not line.startswith('@') and dic_ID_correct.has_key(line.strip().split()[0]):
                print>>output_R2, line.strip()
    output_R2.close()

    return 'R1_forcorrect.{}.sam'.format(suffix), 'R2_forcorrect.{}.sam'.format(suffix)


# --------
class Unmatch:
    '''
    Use to correct the unmatch the Read1 and Read2 pairs:
        Read1 and Read2 are mapped to the same chr and distance of mapped <=5bp

    Generate fastq files for Read1 and Read2 containing corrected sequences
    Name_unmatchcorrect.R1.fastq, Name_unmatchcorrect.R2.fastq
    '''
    
    def __init__(self, samfile_R1, samfile_R2, name, outputdir):
        # both sam files do not have header
        self.samfile_R1 = samfile_R1 # R1 is mapped using bismark
        self.samfile_R2 = samfile_R2 # R2 is mapped using bowtie2
        
        # Run functions to correct the unmatch the Read1 and Read2 pairs
        self.CorrectReadPairs(name, outputdir)

    def ConvertSeq(self, seq, cigar):
        '''
        ConvertSeq for comparison between Read1 and Read2
        Here Skip the Insertion base in return convertedSeq, because insertion base should be removed from mapping most likely
        This is different from AddXMtag.py conversion
        e.g.
        seq =  'TATAAAGAGTTTAGTTTTATATGTTTTTTTTTTTTTTATTTTTTGTATTTATTATGAGTATAAGGTGATTTTCATTTTTAATTATTATTTGATGTTTAT'
        cigar = '23M1I75M'
        return 'TATAAAGAGTTTAGTTTTATATGTTTTTTTTTTTTTATTTTTTGTATTTATTATGAGTATAAGGTGATTTTCATTTTTAATTATTATTTGATGTTTAT'
        '''

        # trim the suspending sequence at both ends, and remove the corresponding cigar
        if 'S' in cigar[0]: # 5'end suspension
            seq = seq[int(re.findall('(\d+)S', cigar[0])[0]):]
            cigar.remove(cigar[0])
        if 'S' in cigar[-1]: # 3'end suspension
            seq = seq[:-1*int(re.findall('(\d+)S', cigar[-1])[0])]
            cigar.remove(cigar[-1])
        
        ls = []
        index = 0
        for item in cigar:
            number = int(re.findall('(\d+)[MID]', item)[0])
            
            if 'D' in item:
                ls.extend(['^'] * number)
            else:
                for i in range(index, index+number):
                    if 'M' in item:
                        ls.append(seq[i])
                index = index + number
        return ''.join(ls)

    def AlignPair(self, R1, R2, start1, start2, strand):
        '''
        Use to align the Read Pair if 5' mapping (Start for forward and End for reverse mapping) is not aligned
        e.g. R1 start=171480982, R2 start= 171480983
        R1: TTATGAAGTTTACTTTGATAGAGATTTGGTGCTATATATAGAGTAGAAATTTAATAGAGTTTTTAGTTGTTATGTTTATTAGATGATATGT
        R2:  CATGAAGTCTACTTTGACAGAGACTTGGTGCTACACATAAAGCAGAAACCCAACAGAGTCTTCAGCTGCTATGTTTACCAGATGATATGT

        return two sequences aligned at 5'end
        R1:  TATGAAGTTTACTTTGATAGAGATTTGGTGCTATATATAGAGTAGAAATTTAATAGAGTTTTTAGTTGTTATGTTTATTAGATGATATGT
        R2:  CATGAAGTCTACTTTGACAGAGACTTGGTGCTACACATAAAGCAGAAACCCAACAGAGTCTTCAGCTGCTATGTTTACCAGATGATATGT
        '''
        distance = start1 - start2
        if strand == '+':
            if distance < 0:
                R1 = R1[abs(distance):]
            else:
                R2 = R2[abs(distance):]
        else:
            if distance > 0:
                R1 = R1[abs(distance):]
            else:
                R2 = R2[abs(distance):]
        return R1, R2

    def ComparePair(self, R1, R2, qual1, qual2):
        '''
        Compare the ConvertedSEQ of R1 and R2 (not including Insertion base) and modify following:
            If a mismatch of R1/R2 is a deletion in one of the read, use the base on the other read and the corresponding qual
            If a mismatch of R1/R2 is a deletion in both reads, pass
        
        Meanwhile adjust the quality score accordingly
        Since Qual does not use symbol ^ as quality score label

        Also shorten the Read Pair to the shortest one
        '''
        ls_R1, ls_qual_R1 = [], []
        ls_R2, ls_qual_R2 = [], []

        for i in range(min(len(R1), len(R2))): # also shorten the Read Pair to the shortest one
            if R1[i] !='^': # R1 base != deletion
                ls_R1.append(R1[i])
                ls_qual_R1.append(qual1[i])
                if R2[i] !='^': # Both R1, R2 base != deletion
                    ls_R2.append(R2[i])
                    ls_qual_R2.append(qual2[i])
                else: # R1 base != deletion but R2 base == deletion -> set R2 base using R1 base
                    ls_R2.append(R1[i])
                    ls_qual_R2.append(qual1[i])
            else: # R1 base == deletion
                if R2[i] !='^': # R1 base == deletion but R2 base != deletion -> set R1 base using R2 base
                    ls_R2.append(R2[i])
                    ls_qual_R2.append(qual2[i])
                    ls_R1.append(R2[i])
                    ls_qual_R1.append(qual2[i])
                else:
                    pass
        return ''.join(ls_R1), ''.join(ls_R2), ''.join(ls_qual_R1), ''.join(ls_qual_R2)

    
    def PairSeq(self, R1_sam, R2_sam):
        '''
        correct the seq of R1 sam entry and R2 sam entry to match R1 and R2
        R1 sam entry and R2 sam entry have the same ID
        '''
        gtf_start_R1 = int(R1_sam.strip().split()[3])
        gtf_start_R2 = int(R2_sam.strip().split()[3])
    
        col10_R1, cigar_R1, col11_R1 = R1_sam.strip().split()[9], re.findall('\d+[MDIS]', R1_sam.strip().split()[5]), R1_sam.strip().split()[10]
        col10_R2, cigar_R2, col11_R2 = R2_sam.strip().split()[9], re.findall('\d+[MDIS]', R2_sam.strip().split()[5]), R2_sam.strip().split()[10]

        gtf_end_R1 = sum([int(re.findall('(\d+?)[MD]', item)[0]) for item in cigar_R1 if 'I' not in item and 'S' not in item]) + gtf_start_R1 - 1
        gtf_end_R2 = sum([int(re.findall('(\d+?)[MD]', item)[0]) for item in cigar_R2 if 'I' not in item and 'S' not in item]) + gtf_start_R2 - 1

        # Use RC of SEQ, reverse cigar and reverse QUAL for reverse mapping
        # So for reverse mapping, the ConvertSEQ has the original direction as the input fastq
        # Do not need to convert back for fastq
        if Flag16(R1_sam.strip().split()[1]):
            col10_R1, col10_R2 = RC_seq(col10_R1), RC_seq(col10_R2)
            cigar_R1, cigar_R2 = cigar_R1[::-1], cigar_R2[::-1]
            col11_R1, col11_R2 = col11_R1[::-1], col11_R2[::-1] # For Flag 16, reverse the QUAL to match RC of col10
            start1, start2 = gtf_end_R1, gtf_end_R2
            strand = '-'
        else:
            start1, start2 = gtf_start_R1, gtf_start_R2
            strand = '+'

        ConvertSEQ_R1 = self.ConvertSeq(col10_R1, cigar_R1)
        ConvertQUAL_R1 = self.ConvertSeq(col11_R1, cigar_R1)
        ConvertSEQ_R2 = self.ConvertSeq(col10_R2, cigar_R2)
        ConvertQUAL_R2 = self.ConvertSeq(col11_R2, cigar_R2)

        # Align the 5'end of mapping if 5'end of read not starting at the same position
        if start1 != start2:
            ConvertSEQ_R1, ConvertSEQ_R2 = self.AlignPair(ConvertSEQ_R1, ConvertSEQ_R2, start1, start2, strand)
            ConvertQUAL_R1, ConvertQUAL_R2 = self.AlignPair(ConvertQUAL_R1, ConvertQUAL_R2, start1, start2, strand)

        correctedcol10_R1, correctedcol10_R2, correctedqual_R1, correctedqual_R2 = self.ComparePair(ConvertSEQ_R1, ConvertSEQ_R2, ConvertQUAL_R1, ConvertQUAL_R2)

        return correctedcol10_R1, correctedqual_R1, correctedcol10_R2, correctedqual_R2
    
    def CorrectReadPairs(self, name, outputdir):
        fastq_R1 = open(os.path.join(outputdir, '{}_unmatchcorrect.R1.fastq'.format(name)), 'w')
        fastq_R2 = open(os.path.join(outputdir, '{}_unmatchcorrect.R2.fastq'.format(name)), 'w')
        with open(self.samfile_R1) as f1, open(self.samfile_R2) as f2:
            R1_sam = f1.readline().strip()
            R2_sam = f2.readline().strip()
            while R1_sam:
                ReadID = R2_sam.split()[0]
                # bismark, e.g. A00336:A00336:HF2VMDRXY:1:1101:10113:20102_1:N:0:CTCGAACAGGAATTGC
                # fastq: ReadID   1:N:0:CAACTCCAACGGATTC, Read ID 2:N:0:CAACTCCAACGGATTC
                sequenceIndex =  ':'.join(R1_sam.split()[0].split('_')[-1].split(':')[1:]) # N:0:CTCGAACAGGAATTGC
                correctedcol10_R1, correctedqual_R1, correctedcol10_R2, correctedqual_R2 = self.PairSeq(R1_sam, R2_sam)
                
                print>>fastq_R1, '@{} 1:{}'.format(ReadID, sequenceIndex)
                print>>fastq_R1, correctedcol10_R1
                print>>fastq_R1, '+'
                print>>fastq_R1, correctedqual_R1

                print>>fastq_R2, '@{} 2:{}'.format(ReadID, sequenceIndex)
                print>>fastq_R2, correctedcol10_R2
                print>>fastq_R2, '+'
                print>>fastq_R2, correctedqual_R2

                R1_sam = f1.readline().strip()
                R2_sam = f2.readline().strip()

        fastq_R1.close()
        fastq_R2.close()
        return fastq_R1, fastq_R2


##-------mainbody
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--Read1', help='sorted input Read1 bismark sam file', dest='Read1_bam')
    parser.add_argument('--Read2', help='sorted input Read2 bowtie sam file', dest='Read2_bam')
    parser.add_argument('--name', help='name for output fastq file', dest='name')
    parser.add_argument('--dir', help='dir path for output file', dest='output_dir', required=False)

    parser.add_argument('--path_to_samtools', help='specify a path to samtools', dest='pathSamtools', default='samtools')

    args = parser.parse_args()

    localtime = time.asctime(time.localtime())
    suffix = localtime.split()[-2].split(':')
    suffix.append(localtime.split()[-1]) 
    suffix= ''.join(suffix)

    # Extract mappable reads and sort by name
    R1_sort_sam, R2_sort_sam = SortBamFile(args.Read1_bam, args.Read2_bam, suffix, args.pathSamtools)
    
    # Extract R1, R2 read pairs for correction
    R1_forcorrect_sam, R2_forcorrect_sam = ExtractReadsForCorrection(R1_sort_sam, R2_sort_sam, suffix)
    
    # Correct the Read pairs
    if args.output_dir:
        outputdir = args.output_dir
    else:
        outputdir = os.getcwd()
    Unmatch(R1_forcorrect_sam, R2_forcorrect_sam, args.name, outputdir)

    ls = [R1_sort_sam, R2_sort_sam, R1_forcorrect_sam, R2_forcorrect_sam]
    for tempfile in ls:
        os.remove(tempfile)
    print "Done."
           




        