# -*- coding: utf-8 -*-
# python2.7
'''
Made on Nov 13, 2020
Modifications on Nov 16, 2020:
    Add function to take and detect compressed Read1 and Read2 files
Modifications on Aug 20, 2022:
    change check_call to shell=False
    change default output file name to name.DeconvolutedRead.fq
    change CheckFormat
'''
'''
Requirement:
Based on python2.7, need pathos.

Purpose:
For Accurateseq paired end reads, compare positions in Read1 and Read2, deconvolution including:
(1) convert Read1 bases to C in deconvoution_R1 if R1=T and R2=C, indicating a Non-methy C conversion during bisulfite
(2) save cytosine methylation status in a report
(3) generate Deconvolution_R1_calibration.fq:
    Calibrate the base quality of R1 bases that are mismatching the R2 bases. 
    The calibrated score are saved in output fastq Deconvolution_Calibration_R1.fq
    This is different from DeconvolutionCalibration_v1.py: v1 generates two fq files having calibrated score for both R1 base and R2 base.

Logic:
    Term:
    R1, R2 mismatch: R1_base != R2_base, excluding the scenario R1=T and R2=C, indicating a Non-methy C conversion during bisulfite
    if len(R1)!=len(R2), the read length is the min(len(R1), len(R2))
    
    Criteria for report read in Deconvolution read:
        total R1, R2 mismatch < 10% of read length (100 base sequencing having 10 mismatches)
        discard read with mismatch >= 10% of read length
    I do not filter read with more than 5mC, which is different from Yi's script.
    
    Deconvolution logic:
    I do not calibrate R1=T, R2=C scenario, e.g. with the following probability:
        Cycle=2 R1=T    R2=C    REF=T   0.5651429445326833      %
        Cycle=2 R1=T    R2=C    REF=C   0.42496460535809316     $

    Deconvolution Read by comparing R1 with R2:
    if R1!=N and R2!=N:
        if R1==R2:
            DeconvolutionBase_R1=R1, DeconvolutionBaseQuality=Quality_R1
            if R1==R2==C, position is marked as a methy C
        else R1!=R2:
            if R1==T and R2==C:
                DeconvolutionBase_R1=C, DeconvolutionBaseQuality=Quality_R1
                position is marked as a Non-methy C
            elif R1==C and R2!=C:
                position is marked as a methy C
                DeconvolutionBase_R1=R1=C, DeconvolutionBaseQuality=Recalibration(REF=R1=C)
            else:
                other R1!=R2 (e.g. R1==G and R2==A)
                DeconvolutionBase_R1=R1, DeconvolutionBaseQuality=Recalibration(REF=R1)
    else:
        if R1==N and R2!=N:
            DeconvolutionBase_R1=R2, DeconvolutionBaseQuality=Quality_R2
            if R2==C, position is marked as a Non-methy C
        else R1!=N and R2==N: Do not change N base in R2
            DeconvolutionBase_R1=R1, DeconvolutionBaseQuality=Quality_R1
            if R1==C, position is marked as a methy C
    Recalibration(REF=R1=X) is calculated by BaseCalibration.py, which is provided in --probability file.

Usage
$python /mnt/home/ettwiller/yan/Bo_script/DeconvolutionCalibration_v2.py \
--Read1 /home/yan/Bo/AccurateSeq/NA12878_AccurateSeq/deconvoluted/AccRep1_R1_val_1.fq \
--Read2 /home/yan/Bo/AccurateSeq/NA12878_AccurateSeq/deconvoluted/AccRep1_R2_val_2.fq \
--name AccRep1 \
--probability /home/yan/Bo/AccurateSeq/200730_AccurateSeq/test/AccRep1.BaseCalibration.probability

--Read1 and Read2: fastq or fastq.gz
    hairpin removed Read1 and Read2, readID should match each other in order.
    if inputs are compressed files, ungzip in the temp dir first.
    The format of input (gz or not) will be automatically detected.
--probability:
    BaseCalibration.probability file generated by BaseCalibration.py, used for Base Quality calibration.
    e.g.
    <Cycle, Read1Base, Read2Base, REF, probability, corresponding quality score>
    Cycle=1 R1=A    R2=T    REF=A   0.674802838049182       &
    Cycle=1 R1=A    R2=T    REF=T   0.3174817528706056      #

--name: name of output files
    output files are saved in dir where the script is run.
    name.DeconvolutedRead.fq: Deconvoluted Read1 with Non-methy C convertd back and sequencing quality calibrated for mismatch bases
    name.Deconvolution.5mC: C or 5mC report

    Note:
    (1) ReadID in output fastq files
    All these three files have modified ReadID, for which the read pair information is removed.
    Original ReadID in Read1: @A00336:A00336:HNKCMDRXX:1:1101:10004:10081 1:N:0:ATTACTCGCCTATCCT
    Original ReadID in Read2: @A00336:A00336:HNKCMDRXX:1:1101:10004:10081 2:N:0:ATTACTCGCCTATCCT
    -> name.Deconvolution.R1.fq ReadID: @A00336:A00336:HNKCMDRXX:1:1101:10004:10081

    (2) name.Deconvolution.5mC
    Here the index starts with 0, corresponding to the first base in Deconvolution_R1: Readseq[0] is the first cycle base.
    C for normal C, M for C resistent to bisulfite conversion.
    readID does not have '@' at the begining
    <readID, methylation information>
    A00336:A00336:HNKCMDRXX:1:1101:10004:10081      C1C8C13C14C29C37C39C45C53C60C61C67C70C71C72C79C80M84C93C94C96
    C1: Readseq[1] is a Non-methy C; M84: Readseq[84] is a methyl C
    Note:
    If there is no C or Methylated C, the methylation if empty, e.g.
    A00336:A00336:HV7F7DRXX:1:1101:1009:7247    
    So the order of reads in name.Deconvolution.5mC is the same as in the name_Deconvolution_Calibration_R1.fq

Note:
(1) Use multiple thread (nodes=ncpu) to speed up, split into sub file containing 50 million reads;
    if pathos.helpers.cpu_count()>8, use nodes=8; else use nodes=pathos.helpers.cpu_count()
(2) tempory dir is named based on current time
'''

try:
    import argparse
    import os, sys
    from subprocess import check_call
    import re
    from pathos.pools import ProcessPool
    import pathos
    script_dir = os.path.dirname(__file__)
    sys.path.append(script_dir)
    from getpath import GetFilePath, CreatePrefix, CheckFormat
except:
    print "Import module error."
    quit()


class BaseCalibrationDict:
    '''
    Use to Extract the probability from BaseCalibration.probability,
    and return a dictionary {(cycle, R1, R2, REF): corresponding quality score}, {('1', A, T, A): &}
    '''
    def __init__(self, input):
        self.input = os.path.abspath(input) # BaseCalibration.probability file
        self.dic_probability = self.ReadReport()

    def ReadReport(self):
        '''
        Use to make a dictionary from probability input
        <Cycle, Read1Base, Read2Base, REF, probability, corresponding quality score>
        Cycle=1 R1=A    R2=T    REF=A   0.674802838049182       &
        
        return {(cycle, R1, R2, REF): corresponding quality score}, {('1', A, T, A): &}
        if probability=NA (e.g. cycle=100), do not include in the dictionary
        '''
        with open(self.input) as f:
            key = [(line.strip().split('\t')[0].split('=')[1], line.strip().split('\t')[1].split('=')[1], line.strip().split('\t')[2].split('=')[1], line.strip().split('\t')[3].split('=')[1]) \
            for line in f if not line.strip().endswith('NA')]
            f.seek(0)
            value = [line.strip().split('\t')[-1] for line in f if not line.strip().endswith('NA')]
            dic_probability = dict(zip(key, value))
            
        return dic_probability
        # {(cycle, R1, R2, REF): corresponding quality score}, {('1', A, T, A): &}
        
class DeconvoluteCalibration:
    '''
    Use to deconvolute by comparing Read1 with Read2,
    calibrate the quality for R1!=R2 base for both R1 and R2,
    extract methylation information
    '''
    def __init__(self, Read1, Read2, dic_probability, name):
        self.Read1 = os.path.abspath(Read1) # hairpin adapter removed Read1
        self.Read2 = os.path.abspath(Read2) # hairpin adapter removed Read2
        self.name = name # name of the output Deconvlution files and methylation file: Deconvolution_R1.name, Deconvolution_R2.name, methylation.name
        self.dic_probability = dic_probability # {(cycle, R1, R2, REF): corresponding quality score}, returned as BaseCalibrationDict.dic_probability
        
        self.count = self.DeconvoluteFastq() # total_input, total_match, total_unmatched, total_more_than_5mC, total_base_matching_length, total_base_mismatch

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
                        ls_R1_qual.append(self.dic_probability.get((str(i), 'C', Read2[i], 'C'), Read1_qual[i])) # {(cycle, R1, R2, REF)
                        ls_C.append('M{}'.format(i))
                        count_mismatch += 1
                    else: # other R1!=R2
                        ls_R1.append(Read1[i])
                        ls_R1_qual.append(self.dic_probability.get((str(i), Read1[i], Read2[i], Read1[i]), Read1_qual[i])) # {(cycle, R1, R2, REF)
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
        output_R1 = open('Deconvolution_R1.{}'.format(self.name), 'w')
        output_methylation = open('methylation.{}'.format(self.name), 'w')
        
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

        return total_input, total_match, total_unmatched, total_more_than_5mC, total_base_matching_length, total_base_mismatch

def main(Read1, Read2, dic_probability, name):
    '''
    Split the large Read1 and Read2 file into file with 100 million reads
    Run class DeconvoluteCalibration for each using multiple threads
    '''
    prefix = CreatePrefix()

    Read1 = GetFilePath(Read1)
    Read2 = GetFilePath(Read2)
    CWD = os.getcwd()
 
    # Using multiple thread
    os.mkdir('Deconvolution{}'.format(prefix)) # Deconvolution172545
    os.chdir('Deconvolution{}'.format(prefix))

    # Check input files are compressed or not
    if CheckFormat(Read1):
        command = 'gunzip -c {} > {}'.format(Read1, 'Read1.fq') # uncompress the gz -> Read1.fastq
        check_call(command, shell=True) # use shell=True for gunzip
        Read1 = 'Read1.fq'
    if CheckFormat(Read2):
        command = 'gunzip -c {} > {}'.format(Read2, 'Read2.fq') # uncompress the gz -> Read1.fastq
        check_call(command, shell=True)
        Read2 = 'Read2.fq'

    # split into 50 million of reads for each thread
    command = ['split', '-l', '200000000', '-d', '-a', '2', Read1, 'R1_'] # 'split -l 200000000 -d -a 2 {} R1_'.format(Read1) # R1_00, R1_01, R1_02 ...
    check_call(command, shell=False)
    command = ['split', '-l', '200000000', '-d', '-a', '2', Read2, 'R2_']
    check_call(command, shell=False)

    ls_input_R1 = []
    ls_input_R2 = []
    for item in os.listdir(os.getcwd()):
        if 'R1_' in item:
            ls_input_R1.append(item)
        elif 'R2_' in item:
            ls_input_R2.append(item)

    if int(pathos.helpers.cpu_count()) > 8:
        n = 8
    else:
        n = int(pathos.helpers.cpu_count())
    print "Use n={} nodes for pathos.".format(n)
    pool = ProcessPool(nodes=n) # nodes=4
    ls_output_name = [item.replace('R1_', '') for item in sorted(ls_input_R1)] # [00, 01, 02 ...]
    ls_return = pool.map(DeconvoluteCalibration, sorted(ls_input_R1), sorted(ls_input_R2), [dic_probability]*len(ls_input_R1), ls_output_name)
    # output files: Deconvolution_R1.00, methylation.00
    # ls_return: a list of call Deconvolute object, each has property count:
        # total_input, total_match, total_unmatched, total_more_than_5mC, total_base_matching_length, total_base_mismatch
   
    pool.close()
    pool.join()

    # combine all the output files 
    with open(os.path.join(CWD, '{}.DeconvolutedRead.fq'.format(name)), 'w') as output:
        for input in ['Deconvolution_R1.{}'.format(x) for x in ls_output_name]: # Deconvolution_R1 for split files
            with open(input) as f:
                for line in f:
                    print>>output, line.strip()
    with open(os.path.join(CWD, '{}.Deconvolution.5mC'.format(name)), 'w') as output:
        for input in ['methylation.{}'.format(x) for x in ls_output_name]: # methylation report for split files
            with open(input) as f:
                for line in f:
                    print>>output, line.strip()

    total_input = sum([item.count[0] for item in ls_return])
    total_match = sum([item.count[1] for item in ls_return])
    total_unmatched = sum([item.count[2] for item in ls_return])
    total_more_than_5mC = sum([item.count[3] for item in ls_return])
    total_base_matching_length = sum([item.count[4] for item in ls_return])
    total_base_mismatch = sum([item.count[5] for item in ls_return])

    print "Deconvolution of Read1 with BaseQuality Calibration."
    print "Read Pairs in input: {}".format(total_input)
    print "Matched Deconvoluted Read Pairs: {}".format(total_match)
    print "Unmatched Read Pairs: {}".format(total_unmatched)
    print "Read having more than 5 5mC: {}".format(total_more_than_5mC)
    print "Deconvoluted bases: {}".format(total_base_matching_length)
    print "Bases of R1!=R2: {}".format(total_base_mismatch)

    print "The output files are:"
    print os.path.join(CWD, '{}.DeconvolutedRead.fq'.format(name))
    print os.path.join(CWD, '{}.Deconvolution.5mC'.format(name))

    os.chdir(CWD)
    command = ['rm', '-r', 'Deconvolution{}'.format(prefix)]
    check_call(command, shell=False)
    return 1
      
##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='hairpin removed Read1 file', dest='Read1')
    parser.add_argument('--Read2', help='hairpin removed Read2 file', dest='Read2')
    parser.add_argument('--probability', help='probability file generated by BaseCalibration.py', dest='probability_file')
    parser.add_argument('--name', help='name suffix used for output', dest='name')
    args = parser.parse_args()

    dic_probability = BaseCalibrationDict(args.probability_file).dic_probability
    main(args.Read1, args.Read2, dic_probability, args.name)
