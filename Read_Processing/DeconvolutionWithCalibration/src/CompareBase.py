# -*- coding: utf-8 -*-
# python2.7
'''
Made on Oct 16, 2020
Modifications on Oct 26, 2020:
    Add (if Read[index] != '^') to skip the deletion base, previously cycle could be more than 100 due to mishandling of cigar D
    Skip the RemoveKnownVariant step and use all the bases for analysis if vcf file is not provided.
Modifications on Oct 29, 2020:
    Change starting cycle to 0 to match the python index
Modifications on Jan 12, 2021:
    Change to couple the convertedSEQ and convertedSEQ_Index using a dictionary to better calculate the base cycle
Modifications on Aug 3, 2021:
    Add BAQ and MAPQ cutoff option, default 0
    Change to work for both single or pair end read ID in SaveRef()
Modifications on Aug 20, 2022:
    change check_call to shell=False
    add option to specify bedtools and samtools
'''
'''
Requirement:
Based on python2.7, need samtools, bedtools executable

Purpose:
Use to compare the base in read with the corresponding genome sequence.
The output could be used for base quality calibration.

Logic:
    I use only the forward mapping for error rate estimation since the error rate between forward and reverse mapping is almost the same.
    The bam file (downsample decided by percent) is converted to bed and get corresponding genome sequence using getfasta;
    compare each base in read in the fastq with the corresponding genome base,
    if read in fastq does not exist in bam (downsample bam), skip this read;
    At the end remove the know SNP sites if vcf file(s) is provided, also apply BAQ and MAPQ cutoff.
    Note:
    In principle can extract the read/BAQ information from sam file (converted from bam) therefore do not need the fastq file.
    Since use dict to save the bam/getfasta information, reads in fastq and bam do not need to have the same order.

Explanation:    
    $samtools view -h AccRep1_R2_default.uni.nodup.bam -s 0.0005 -b > test.bam
    $bedtools bamtobed -cigar -i test.bam > test.bed
    $bedtools getfasta -fi /media/neb/backupHD/GCA_000001405.15_GRCh38_no_alt_spike.fa -bed test.bed -s -bedOut > test.fa

    ^ Reads mapped to reverse strand:
    test.bam:
    A00336:A00336:HNKCMDRXX:2:2217:23176:17174      16      chr1    823780  6       90M     *       0       0       AATGAGGATATAAAAATGACATGAGTAAAGCCTGAAAATCATATTATTAGAATGGAGGACATAGTGACATAAAACCTCAGTGACTGTAGA      FFFFFFFFFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      AS:i:-3 XS:i:-8 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:24A65      YT:Z:UU
    test.fa: <chr, strand, end, ID, MAPQ, strand, cigar, genome seq>
    chr1	823779	823869	A00336:A00336:HNKCMDRXX:2:2217:23176:17174	6	-	90M	TCTACAGTCACTGAGGTTTTATGTCACTATGTCCTCCATTCTAATAATATGATTTTCAGGCTTTATTCATGTCATTTTTATATCCTCATT
    Bottom Genome: TCTACAGTCACTGAGGTTTTATGTCACTATGTCCTCCATTCTAATAATATGATTTTCAGGCTTTATTCATGTCATTTTTATATCCTCATT
    Fastq seq:     TCTACAGTCACTGAGGTTTTATGTCACTATGTCCTCCATTCTAATAATATGATTTTCAGGCTTTACTCATGTCATTTTTATATCCTCATT
    Here ReadBase=C,REF=T (REF is the base on bottom genome), so this misreading of T(REF) to C(ReadBase) results in REF=A,ALT=G for variant calling

    ^ Reads mapped to forward strand:
    test.bam:
    A00336:A00336:HNKCMDRXX:2:1223:2944:23938       0       chr1    1281698 42      100M    *       0       0       TCTGACTGAGGCCCTGCAGGGCATACGGGTCATGGAAGGGGTGCTCTGCCCCCGGCCAGGACTCCCCCTCCTACAGGAAGCCGGGGGCCTTGCTCCTGCA    FFFFFFFFFFFFFFFFF:FFFFFFFFFFF:F:FFFFFFFFFFFF,FFFFFFFFFFFFFFF:FFFFFFF:FFF,FFFFF:FFFFFFFFFFFFFFFFFFFFF    AS:i:-4 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:60G39      YT:Z:UU
    test.fa:
    chr1	1281697	1281797	A00336:A00336:HNKCMDRXX:2:1223:2944:23938	42	+	100M	TCTGACTGAGGCCCTGCAGGGCATACGGGTCATGGAAGGGGTGCTCTGCCCCCGGCCAGGGCTCCCCCTCCTACAGGAAGCCGGGGGCCTTGCTCCTGCA
    Top Genome: TCTGACTGAGGCCCTGCAGGGCATACGGGTCATGGAAGGGGTGCTCTGCCCCCGGCCAGGGCTCCCCCTCCTACAGGAAGCCGGGGGCCTTGCTCCTGCA
    Fastq seq:  TCTGACTGAGGCCCTGCAGGGCATACGGGTCATGGAAGGGGTGCTCTGCCCCCGGCCAGGACTCCCCCTCCTACAGGAAGCCGGGGGCCTTGCTCCTGCA
    Here ReadBase=A,REF=G (REF is the base on top genome), so this misreading of G(REF) to A(ReadBase) results in REF=G,ALT=A for variant calling

    Conclusion: 
    forward or reverse mapping has the same kind of error rate: forward mapping compared to top strand ssDNA while reverse mapping compared to bottom strand ssDNA,
    note this is different from variant calling in which mapping strand matters.
    See ConvertSeq and ExtractRead steps for details about how to compare Read base to reference with consideration of deletion and insertion.

Usage:
$python CompareBase.py 
--bam TestSeq.bam
--vcf /home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf 
--fastq TestSeq_R1_val_1.fq 
--output TestSeq_deconvoluted_R1.compareBase.txt
--percent DefaultValue
--reference DefaultValue
--BAQ DefaultValue
--MAPQ DefaultVale

--bam:
    bam or sam file of Read1 (or Read2, but can not contain both Read1 and Read2) corresponding to the --fastq file
    This script will remove unmapped reads (-F 4), so do not need to apply -F 4 on the bam file before this step.
--fastq:
    fastq or fastq.gz of original reads removing adapter sequence, which is the file used for mapping.
--output: 
    ReadID, mapq, cycle, chr, 1-coordination position, Reference_base(REF), Read_base, Read_baseQuality(score not symbol)
    A00336:A00336:HNKCMDRXX:1:1104:29044:29684      42      0       chr4    148667385       G       G       37
    Note:
    cycle starting from 0 (0-99 for 100bp sequencing) to match the Readseq index.
    since only use forward mapping, the Reference_base(REF) represents the top genome base.
--reference: Default /mnt/home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/GCA_000001405.15_GRCh38_no_alt_spike.fa
    Reference genome used for mapping to generate the bam file

--percent: float, Default 1
    percent for downsample using: samtools view -h -s percent
    if percent=1, skip the samtools downsample step, use all the reads in bam file for analysis 

--vcf: Default None 
    One vcf file containing known polymorphic sites; the known polymorphisms are not saved in output.
    e.g. /mnt/home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf
    if not provided, use all the bases for analysis
    Several variant files can be merged into one using picard.jar MergeVcfs function.

--BAQ: int [0, 40], Dafault 0
    Read Base quality cutoff, ASCII-33=BAQ, only save Base with BAQ>=cutoff in output, filter based on col8
--MAPQ: int, Dafault 0
    Read mapping quality cutoff, only save Base from reads with MAPQ>=cutoff in output, , filter based on col2

--path_to_bedtools, --path_to_samtools: Optional
    use this option to specify a path to bedtools/samtools executable, e.g.
    /usr/bin/bedtools
    Else it is assumed that bedtools/samtools is executable in the PATH.

Note:
(1) Need to assign enough memory on cluster for the ExtractionReference step
    memory useage: need to create a dict saving {ReadID: (chr, start (1-coordination), cigar, mapq, reference seq)} for all the reads
(2) Only choose the reads mapped to forward strand (-F 16) for analysis for convenience,
    samtools view -F 4 and samtools view -F 16 are applied before bamtobed step.
    As shown above, the error rate between ReadBase and REF (bottom genome for reverse mapping) is the same for forward and reverse mapping.
    Therefore I use only the forward mapping for error rate estimation.
    In addition, when recalibrate the quality by comparing Read1 to Read2 in Accurateseq, I only need to consider the error type but not mapping strandness.
(3) Exclude reads with suspendsion (cigar S) from analysis.
(4) cycle starts at 0, so the cycles for 100bp sequencing is (0, 1 ...99), which matches the python index of read. Not all the reads have full cycles.
(5) tempory files are named based on time, could be called by different runs.
(6) size of tempory files
    5.5G	AccRep3_lane1_Deconvolution_R1.fq; 7.1G	AccRep3_lane1_deconvoluted_R1.sam -> 60G prefix.temp.output
'''


import re
import os, sys
import time
from subprocess import check_call
import argparse
import gzip

script_dir = os.path.dirname(__file__)
sys.path.append(script_dir)
from getpath import CheckFormat


class ExtractionReference:
    '''
    Downsample bam (percent=0.1) and extract forward strand mapping
    Convert to bed and fasta to get reference sequencing
    e.g. Use AccRep2 forward strand mapping as examples:
    samtools view -h -s 0.1 /home/yan/Bo/AccurateSeq/200730_AccurateSeq/deconvoluted/bowtie2/AccRep2_deconvoluted_default.uni.nodup.bam | samtools view -h -F 16 - -b > test_forward.bam
    bedtools bamtobed -cigar -i test_forward.bam > test_forward.bed
    bedtools getfasta -fi /home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/GCA_000001405.15_GRCh38_no_alt_spike.fa -bed test_forward.bed -s -bedOut > test_forward.fa
    test_forward.fa e.g. chr, start, end, ReadID, MAPQ, strand, cigar, reference
    
    Extract information from bedtofasta fasta file and save in a dictionary:
        {ReadID: (chr, start (1-coordination), cigar, mapq, reference seq)}
    '''
    def __init__(self, bam, percent, REF, pathBedtools, pathSamtools):
        localtime = time.asctime(time.localtime())
        self.prefix = ''.join(localtime.split()[-2].split(':')) # '151542'
        self.bam = os.path.abspath(bam) # mapping bam file corresponding to the Read.fastq
        self.percent = percent
        self.REF = REF # reference geome fa file used for mapping and bedtools getfasta

        self.pathBedtools = pathBedtools # /usr/bin/bedtools 
        self.pathSamtools = pathSamtools
        
    def GetFa(self): # percent: samtools view -s percent, default 0.1
        '''
        make sure have enough memory for samtools step
        '''
        with open('CompareBase_{}.sh'.format(self.prefix), 'w') as output:
            print>>output, '#!/bin/sh'
            # print>>output, 'REF=\'/home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/GCA_000001405.15_GRCh38_no_alt_spike.fa\''
            print>>output, 'REF={}'.format(self.REF)
            if self.percent != 1:
                command = '{} view -h -s {} {} | {} view -h -F 4 - -b > {}_map.bam'.format(self.pathSamtools, \
                    self.percent, self.bam, self.pathSamtools, self.prefix)
                print>>output, command
            else:
                command = '{} view -h -F 4 {} -b > {}_map.bam'.format(self.pathSamtools, self.bam, self.prefix)
                print>>output, command
            command = '{} view -h -F 16 {}_map.bam -b > {}_forward.bam'.format(self.pathSamtools, self.prefix, self.prefix) # prefix_forward.bam
            print>>output, command
            command = '{} bamtobed -cigar -i {}_forward.bam > {}_forward.bed'.format(self.pathBedtools, self.prefix, self.prefix) # prefix_forward.bed
            print>>output, command
            command = '{} getfasta -fi $REF -bed {}_forward.bed -s -bedOut > {}_forward.fa'.format(self.pathBedtools, \
                self.prefix, self.prefix) # prefix_forward.fa
            print>>output, command
        try:
            check_call(['sh', 'CompareBase_{}.sh'.format(self.prefix)], shell=False)
        except:
            print 'Error in the convertion of bam to get fasa reference.'
            quit()

        # remove temp files 
        ls = ['{}_map.bam'.format(self.prefix), '{}_forward.bam'.format(self.prefix), '{}_forward.bed'.format(self.prefix), 'CompareBase_{}.sh'.format(self.prefix)]
        for item in ls:
            os.remove(item)
        return '{}_forward.fa'.format(self.prefix)

    def SaveRef(self):
        '''
        Create a dict saving the ReadID, cigar and corresponding reference sequence
        pair end, getfasta return ID A00336:A00336:HHKJFDRXY:1:2103:4689:14653/1, ID in fastq: @A00336:A00336:HHKJFDRXY:1:2101:10004:10207 1:N:0:GCCTTAACTAGGAGCT
        single end, getfasta return ID A00336:A00336:HHKJFDRXY:1:2103:4689:14653, ID in fastq: @A00336:A00336:HHKJFDRXY:1:2101:10004:10207 1:N:0:GCCTTAACTAGGAGCT
        '''
        with open(self.reference) as f:
            key = [line.strip().split()[3].split('/')[0] for line in f] # ID has the same fortmat as in fastq
            f.seek(0)
            value = [(line.strip().split()[0], int(line.strip().split()[1])+1, line.strip().split()[6], line.strip().split()[4], line.strip().split()[7]) for line in f]
            # {ReadID: (chr, start (1-coordination), cigar, mapq, reference seq)}
            reference_dic = dict(zip(key, value))
        return reference_dic
    
    def GetDict(self):
        '''
        Run two functions to get reference dictionary
        '''
        self.reference = self.GetFa() # samtools view -s percent
        self.reference_dic = self.SaveRef() # {ReadID: (chr, start (1-coordination), cigar, mapq, reference seq)}
        return self.reference_dic


def ConvertSeq(seq, cigar):
    '''
    Use to return read sequence matching the Cigar, compared with reference sequence

    Logic:
    use '^' to represent bases corresponding to D (deletion compared to reference) in cigar.
    use '-' to represent bases corresponding to I (insertion compared to referece) in cigar.

    seq = 'GTGATGCAGCTCTTCGGCCTGGTTAACACCCTTCTGGCCAATGACCCAACATCCCTTCGGAAAANCTC'
    cigar = ['60M', '2I', '2M', '1D', '4M']
    return GTGATGCAGCTCTTCGGCCTGGTTAACACCCTTCTGGCCAATGACCCAACATCCCTTCGG--AA^NCTC

    seq = 'ATAGACTACATCTANATAATTGCCTCTTCAAATAGTCGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT'
    cigar = ['36S', '64M']
    'CGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT'

    seq = 'CCTAACCCTANCCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCNAACCNTAACNCTA'
    cigar = ['41M', '1I', '41M', '1D', '17M']
    'CCTAACCCTANCCCTAACCCTAACCCTAACCCTAACCCTAA-CCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA^CCCNAACCNTAACNCTA'
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
                else: # I in item
                    ls.append('-')
            index = index + number
    return ''.join(ls)

def FindREF_Index(convertedSEQ):
    '''
    return a dict {convertedSEQ_Index: REF_Index}
    '''
    ls_convertedSEQ_Index = [i for i in range(len(convertedSEQ))]
    i = 0
    ls_REF_Index = []
    for item in convertedSEQ:
        if item != '-':
            ls_REF_Index.append(i)
            i += 1
        else:
            ls_REF_Index.append('-')
    dic2_REF_Index = dict(zip(ls_convertedSEQ_Index, ls_REF_Index)) # {convertedSEQ_Index: REF_Index}
    return dic2_REF_Index


def ExtractRead(reference_dic, input_fastq, output_file):
    '''
    Extract read base information and corresponding reference base
    output file:
    ReadID, mapq, cycle, chr, 1-coordination position, Reference_base, Read_base, Read_baseQuality
    '''
    output = open(output_file, 'w')
    total_count = 0 # total number of reads used for analysis
    # Deconvoluted_R1 or R2 fastq file, choose reads listed in reference_dic and do not have suspension (cigar S) for analysis
    if CheckFormat(input_fastq):
        f = gzip.open(input_fastq)
    else:
        f = open(input_fastq)

    line1 = f.readline().strip()
    line2 = f.readline().strip()
    line3 = f.readline().strip()
    line4 = f.readline().strip()
    while line1:
        readID = line1.strip().split()[0].replace('@', '')
        # R2 @A00336:A00336:HV7F7DRXX:1:1101:10004:10019 2:N:0:ATTACTCGCCTATCCT -> A00336:A00336:HV7F7DRXX:1:1101:10004:10019
        # Deconvolution R1 @A00336:A00336:HV7F7DRXX:1:1107:24044:3803 -> A00336:A00336:HV7F7DRXX:1:1107:24044:3803

        # reference_dic has only the forward mapping reads
        if reference_dic.has_key(readID): 
            cigar = re.findall('\d+[MDIS]', reference_dic[readID][2])
            if 'S' not in reference_dic[readID][2]: # cigar does not have suspension
                total_count += 1
                start, reference = reference_dic[readID][1], reference_dic[readID][4] # reference sequencing
                Read = ConvertSeq(line2, cigar) # convert Read
                Read_qual = ConvertSeq(line4, cigar)
                dic2_REF_Index = FindREF_Index(Read) # {convertedSEQ_Index: REF_Index}
                cycle = 0
                for index in range(0, len(Read)):
                    if Read[index] != '^': # skip the deletion base and Insertion base, do not report, but '-' base is included for cycle counting
                        if Read[index] != '-':
                            REF_position = start + dic2_REF_Index[index] # REF position = POS(col4) + REF_Index
                            # ReadID, mapq, cycle, chr, 1-coordination REF position, Reference_base, Read_base, Read_baseQuality
                            ls = [readID, reference_dic[readID][3], str(cycle), reference_dic[readID][0], str(REF_position), reference[dic2_REF_Index[index]], Read[index], str(ord(Read_qual[index])-33)]
                            print>>output, '\t'.join(ls)
                        cycle += 1
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()
        line4 = f.readline().strip()
    output.close()
    f.close()
    print "There are {} number of reads used for Base Calibration analysis.".format(total_count)
    return output_file

def RemoveKnownVariant(vcf, input_file, output_file, BAQ, MAPQ):
    '''
    For positions with Read_base != Reference_base, if it is a known SNP position, remove it
    Apply BAQ and MAPQ cutoff
    input_file is the output of ExtractRead().
    ReadID, mapq, cycle, chr, 1-coordination position, Reference_base, Read_base, Read_baseQuality
    '''
    # remove KNOW variant sites defined in ls_vcf files
    # remove positions if reference is 'N'
    # remove positions below BAQ, MAPQ cutoff
    
    count_total = 0 
    count = 0
    output = open(output_file, 'w')

    if vcf :
        variant_dic = {} # {(chr, 1-coordination position)}
        with open(vcf) as f:
            key = [(line.strip().split()[0],line.strip().split()[1]) for line in f if not line.startswith('#')]
            variant_dic.update(dict.fromkeys(key)) 
        with open(input_file) as f:
            for line in f:
                count_total += 1
                if not variant_dic.has_key((line.strip().split()[3], line.strip().split()[4])) and line.strip().split()[5] != 'N':
                    if int(line.strip().split()[1]) >= MAPQ and int(line.strip().split()[7]) >= int(BAQ): # ord('BAQ Symbol')-33=BAQ (0 to 40)
                        print>>output, line.strip()
                        count += 1
    else:
        with open(input_file) as f:
            for line in f:
                count_total += 1
                if line.strip().split()[5] != 'N':
                    if int(line.strip().split()[1]) >= MAPQ and int(line.strip().split()[7]) >= int(BAQ):
                        print>>output, line.strip()
                        count += 1
    output.close()
    print "There are {} bases in the input fastq file.".format(count_total)
    print "There are {} bases that are not SNP sites and above cutoff BAQ {} and MAPQ {}, saved in the output file.".format(count, BAQ, MAPQ)


def main(dic_args):
    '''
    Provide information and Run function.
    Here I use dic_args to faciliate running this function as an imported module instead of from command line
    dic_args = {'input_bam': , 'percent': ,'REF': , 'input_fastq':, 'vcf':, 'output':, 'BAQ':, 'MAPQ': 'pathBedtools': 'pathSamtools':}
    '''
    print "Compare the Read Base to reference genome."
    print "The bam file used is {}, and downsample to {}.".format(dic_args['input_bam'], dic_args['percent'])
    print "The fastq file used is {}".format(dic_args['input_fastq'])

    temp = ExtractionReference(dic_args['input_bam'], dic_args['percent'], dic_args['REF'], dic_args['pathBedtools'], dic_args['pathSamtools'])
    prefix = temp.prefix
    reference_dic = temp.GetDict()

    temp_output = ExtractRead(reference_dic, dic_args['input_fastq'], '{}.temp.output'.format(prefix)) # temp_output: CompareBase.txt having all the bases.

    # remove the known snp sites and apply BAQ or MAPQ cutoff
    RemoveKnownVariant(dic_args['vcf'], temp_output, dic_args['output'], dic_args['BAQ'], dic_args['MAPQ'])
    
    os.remove(temp_output)
    os.remove('{}_forward.fa'.format(prefix))

##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', help='Read mapping bam file', dest='input_bam')
    parser.add_argument('--fastq', help='Read fastq file', dest='input_fastq')
    parser.add_argument('--output', help='output file', dest='output')
    parser.add_argument('--percent', help='downsample percent', dest='percent', type=float, default=1.0)
    parser.add_argument('--reference', help='reference genome fa file used for mapping', dest='REF', default='/mnt/home/yan/Bo/AccurateSeq/hg38_lambda_XP12_T4_pUC19_combine/GCA_000001405.15_GRCh38_no_alt_spike.fa')
    
    parser.add_argument('--vcf', help='vcf file showing snp', dest='vcf', required=False, default = None)    
    parser.add_argument('--BAQ', help='Read Base quality cutoff', dest='BAQ', type=int, default=0)
    parser.add_argument('--MAPQ', help='Read mapping quality cutoff', dest='MAPQ', type=int, default=0)

    parser.add_argument('--path_to_bedtools', help='specify a path to bedtools', dest='pathBedtools', default='bedtools')
    parser.add_argument('--path_to_samtools', help='specify a path to samtools', dest='pathSamtools', default='samtools')

    args = parser.parse_args()

    main(vars(args))


   