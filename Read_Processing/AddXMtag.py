# -*- coding: utf-8 -*-
# python2.7
'''
# Logs:

Created on March 24, 2020

Modifications on April 23, 2020:
    For entry for which the methylation information is empty in deconvolution report (read containing no C),
    instead of discard the read, report XM tag using '.' for all the positions;
    therefore prefix.noXMtag.sam file only contains the entry for which the ID is not included in deconvolution report.
Modifications on May 7, 2020:
    Add options to split the Read1 sam and report to run multiple threads for a large input file.
Modifications on July 23, 2020:
    Change to use the basename of self.sam to create bamtobed script.
Modifications on Aug 20, 2020:
    For multiple thread, split the sam file without header, otherwise Addtagsam_00 has duplicated header.
    Fix the bug when cigar has 'S'
    Split GetFasta step into multiple steps instead of using stdin to avoid process killed due to memory usage.
Modifications on Nov 17, 2020:
    For split files with --thread option, sort the split file based on name (add sam_list_temp) to ensure that the output has the same order as input.
Modifications on Jan 5, 2020:
    Use bedtools slop -s -r 2 -l 0 to get seq including 2bp downstream. 
    This will be faster than using command to 'split based on strand and add 2bp downstream seperately'
    I compared the two results. They are the same.
Modifications on June 15, 2021:
    Since newseq in Convert2XM() function is not used, do not execuate the conversion for newseq
Modifications on Feb 2 , 2022:
    Add option to change the fa and fa.fai for reference genome and samtools path
Modifications on Aug 20 , 2022:
    Use GetFilePath function to replace os.path.abspath().
    Add --smp to control the nodes
'''
'''
Requirement:
(1) Based on python2.7, need pathos module to run multiprocessing
(2) Bedtools and samtools; command 'bedtools' and 'samtools' should be executable in $PATH, otherwise need to speficy the path in Getfasta step

Purpose:
Use to Generate a XM tag to each entry in the sam file of the Deconvoluted Read mapping;
The determination of methylation status follows the bismark principle based on the methylation report.
The methylation status at each position could be extracted using bismark extractor using output sam file as input.

Usage:
$python AddXMtag.py --input TestSeq_DeconvolutedRead.uni.nodup.sam \
    --report TestSeq.Deconvolution.all.5mC \
    --name TestSeq \
    --reference hg38.fa

--input:
    A sam file that is generated by aligning the Deconvoluted Read to the reference genome.
    AddXMtag.py performs -F 4 to remove the unmapped reads, but not apply other Flag tag filter,
    therefore need to perform -F 256 or -F 2048 before this step when necessary.
    Currently the maximum number of mapping of input is 2970 million if add --thread option to run multiprocessing.

--report:
    The methylation report generated during Deconvolution Step.
    e.g.
    <ReadID, Methylation>
    M03991:M03991:000000000-C7G4V:1:2117:16417:8301	C0C3C10C17C21C23C32C42M46

--name: name for Output files, Default addtag
    e.g. --name TestSeq
    TestSeq.XMtag.sam: saves all the entries with XMtag added;
    TestSeq.noXMtag.sam: saves all the entries that can not be added XMtag: 
    e.g. unmapped reads, reads without information in deconvolution report or reads can not be analyzed for methylation context (at the boundary of chr).
    
    The reads in output are in the same order as in input sam file,
    so if the input sam is sorted, the output sam is also sorted.

--dir: Optional
    Directory to save the output files, not required. 
    If not provided, output files are save in the current working dir.
    
--thread: Optional
    Add --thread to run multiple threads using pathos if input Read1.sam is a large file. 
    The Read1.sam file is split to have 30 million aligments for each using command: 'split -l 30,000,000 -d -a 2 Read1.sam Addtagsam_'
--smp: Optional Int
    number of nodes used for pathos multiple processing.
    If --thread is USED but not define --smp, use all the available CPUs with a limit 8.

--reference: reference.fa needed for bedtools getfasta, which is the reference genome used for alignment.

--pathSamtools, --pathBedtools:
    Use this option to specify a path to the samtools or bedtools executable, 
    e.g. /usr/bin/samtools or /usr/bin/bedtools. 
    Else it is assumed that samtools and bedtools is in the PATH.

Logic:
    Creat XM tag for methylation status corresponding to the position in SEQ;
    the suspending part in SEQ (S in cigar) is not shown in XM tag, so len(XM) = sum of cigar M/I;
    Do not consider the scenario that cigar has hard clip H. 
    See BedTools_README for detailed explainations about comparison between SEQ and corresponding genome sequence.
    See Bismark_README.txt for details about rules of bismark adding XM tags.
    
    output prefix.noXMtag.sam:
        Have discarded reads: 
            Reads do not have methylation information in report, 
            or not aligned (Flag 4) are saved in output, 
            or at the end of chromsome that genome context is not available for all the positions in read,
            or readID is not included in deconvolution report.
        prefix.noXMtag.sam will be empty if no reads discarded.
    
    output prefix.XMtag.sam (not sorted):
        Have Aligned Reads (-F 4) with methylation information (XM tag).

    To report XM tag based on the deconvolution report:
        nonC: position is not a C or Methylated C

        If a position is a C (or Methylated C) in the deconvolution report, and also C on genome:
            Find the methylation context based on the seq in genome, use x,h,z,u or Capital in XM tag following bismark rules.
        If a position is nonC in either deconvolution report or genome:
            Use '.' in XM tag in the corresponding position in SEQ
        If a position is a C on genome but is nonC in deconvolution report:
            Use '.' in XM tag.
        If a position is C in deconvolution report, but nonC on genome (mismatch) or this base is a insertion (I in cigar) not present in genome:
            Use '.' in XM tag.
        
        The C context is determined based on the corresponding genome, in other words based on the mapping strand from 5' to 3' direction:
            For read mapping to reverse strand (Flag 16), 
            e.g. Read in fastq 5'-GcGT-3', SEQ=ACgC (corresponding to top strand 5'->3' direction), bottom genome=5'-GcGT-3', top genome='5-ACgC-3',
            the context is cG depending on the bottom genome .z.., so the XM tag corresponding to SEQ in bam file is reverse of '.z..' -> '..z.' 
        
        So for C at mismatching/SNP site or before mismatching/SNP site, the C context calling could be not accurate.
    
    There is Difference between AddXMtag.py and bismark for calling C coming in front of a deletion (^):
    AddXMtag.py decides the context based on the genome context:
            col10:  	TTTC^^AAATTATTTTGTGATGTGTGTGTTTAATTTATAGAGTTTAATTTTTTTTTTTATAGGGTAGTTTGGAAATAT
            my XM:      .h.X  ...h..h.............z...h..h.h.x.........hh...h....h......x..........h.h
        top genome:   attctcagaaactactttgtgatgtgtgcgttcaactcacagagtttaacctttcttttcatagggcagtttggaaacactc
    bismark labels the C in front of deletion depending on the context of sequence (CAA -> H):
            col10:  	TTTC^^AAATTATTTTGTGATGTGTGTGTTTAATTTATAGAGTTTAATTTTTTTTTTTATAGGGTAGTTTGGAAATAT
        bismark XM:  	.h.H  ...h..h.............z...h..h.h.x.........hh...h....h......x..........h.h
        top genome:   attctcagaaactactttgtgatgtgtgcgttcaactcacagagtttaacctttcttttcatagggcagtttggaaacactc
    
    No Difference between AddXMtag.py and bismark for calling C coming in front of an insertion (-):
    AddXMtag.py decides the context based on the genome context:
            col10:  	TTTC-TTTTTTTTGGGCTTTAGTTTCTTTTTTGGTAAAACGGGGATGGTAATGGGATATTCTTAGGGGGTGTATGA
            my XM:  	...H.....hhx....Hh.x.....Hh.h..........Z................h.hhH.x........h....
        top genome:   gatttc atttccctgggcctcagtttcctctttggtaaaacggggatggtaatgggacaccctcagggggtgcatgagg
    bismark labels the C in front of insertion based on the genome context (CAT -> H):
            col10:  	TTTC-TTTTTTTTGGGCTTTAGTTTCTTTTTTGGTAAAACGGGGATGGTAATGGGATATTCTTAGGGGGTGTATGA
        bismark XM:  	...H.....hhx....Hh.x.....Hh.h..........Z................h.hhH.x........h....
        top genome:   gatttc atttccctgggcctcagtttcctctttggtaaaacggggatggtaatgggacaccctcagggggtgcatgagg

    No Difference between AddXMtag.py and bismark for calling C at or coming in front of mismatch.
    
    Bismark tags needed for bismark extractor:
        XM-tag (methylation call string)
        XR-tag (read conversion state for the alignment)
        XG-tag (genome conversion state for the alignment)

Note:
(1) 
Create a temporary dir AddtagTmp{time} based on current time, so AddXMtag.py could be called by multiple runs at the same time.

(2)
The deconvolution report may contain entries without methylation information (col2 is empty),
which means there is no C in the original fastq read (or no G in Read2).
e.g. 
^ AccTar_deconvoluted_5mC.txt
M03991:M03991:000000000-C7G4V:1:2114:18720:7168
^ AccTar_deconvoluted.fq
@M03991:M03991:000000000-C7G4V:1:2114:18720:7168
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATANA
In this case, report XM tag using '.' for all the positions.

(3) Can only add XM tags to Deconvoluted Read.

(4) Flag tag filter 
Since bowtie2 reports unmapped reads in sam file,
AddXMtag.py performs -F 4 to remove the unmapped reads, but not apply other Flag tag filter,
therefore need to perform -F 256 or -F 2048 before this step when necessary.

(5) Add --thread to speed up for a large sam file, memory and space requirement
Since Split and extraction entries from deconvolution report steps also take time, 
using multiple threads does not linearly increase the speed.

Space for temporary files:
To run multiple thread, a large sam file is split into each having 30 million entries.
The size of a sam file having 30 million entries is 8.4G for 100bp sequencing.
Memory for multiple threads:
Need to assign enough memory when running multiple threads using pathos:
    I run pathos.ProcessPool(nodes=8), which has maxvmem=146G.

Timing:
It needs about 5h to finish for an input sam file having 470 million alignments with node=8.

Currently the maximum number of mapping of input is 2970 million (30 million * 99 split files) if add --thread option to run multiprocessing.

(6) output sam files sorting coordinate
The reads in output files keep the same order as in input sam file,
if the input sam is sorted, the output sam is also sorted.
bismark extractor does not require sorted file.
'''

import re
import numpy as np
from subprocess import check_call
import argparse
import os
import time
import datetime
try:    
    from pathos.pools import ProcessPool
    import pathos
    multipleThreads = True
except:
    print "pathos module error, can not run multiple threads."
    multipleThreads = False

__version__ = '2022.08.20'

##########
def Flag16(Flag):
    '''
    return True if Flag contains Flag 16, Flase if not
    Flag 0 return False
    '''
    if len(str(np.binary_repr(int(Flag))))<5: # Flag < 16
        return False
    else:
        return str(np.binary_repr(int(Flag)))[-5] == '1'

def Flag4(Flag):
    '''
    return True if Flag contains Flag 4, Flase if not
    '''
    if len(str(np.binary_repr(int(Flag))))<3: # Flag < 4
        return False
    else:
        return str(np.binary_repr(int(Flag)))[-3] == '1'

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

##########
class MethylationExtraction:
    def __init__(self, report, sam, fa):
        self.report = report # deconvolution report
        self.sam = sam # deconvoluted Read1 alignment sam file
        self.fa = fa # genome reference
        
        self.contextgenome = {} # 2bp downstream for context genome for each read ID, {read ID: context genome ...}
        self.metylationstatus = {} # methylation status for each ID based on report, {read ID: C4C6C8C10C11C22C24C29C35C37C40C41C43C45C46C48C49C51C52C53C55C56C58}

    def MethylationStatus(self):
        '''
        make a dict saving methlyation status for each ID from report
        
        There are IDs without methylation status in report (there is no C in the original fastq read),
        for these cases use -1 as value, do not use None to distinguish the case that ID does not exist in report.
        '''
        print "Reading Deconvolution methylation report."
        with open(self.report) as f1:
            key = [line.strip().split()[0] for line in f1]
            f1.seek(0)
            value = [line.strip().split()[1] if len(line.strip().split())==2 else -1 for line in f1]
            self.metylationstatus = dict(zip(key, value))
        
        total_noInformation = 0
        for item in self.metylationstatus:
            if self.metylationstatus[item] == -1:
                total_noInformation+=1

        print "There are {} entries having no methylation information in deconvolution report.".format(total_noInformation)
        # print "Memory usage by dict containing methylation status: {}.".format(sys.getsizeof(self.metylationstatus)) 
        return 1
    
    def GetFasta(self):
        '''
        Convert sam to bed, 
        then get fasta of sequence, for which adding 2bp downstream for context determination

        Here need to consider the mapping strandness since methylation centext and downstream sequence depend on strand.
        Therefore I deal with forward and reverse mapping seperately.

        bedtoosl getfasa -s returns the context genome (top genome for Flag 0 and bottom genome for Flag 16) with 2bp downstream
        e.g.
        chr9	69864283	69864365	M03991:M03991:000000000-C7G4V:1:1101:10308:12324_1:N:0:GCCAAT	0	-	45M2D31M	ttctgAAAAAAGcaaacacccacacacacacacacgtatacacacatgcacatacacacacacacacaATTTCCAAAAGCTA
        chr1    16760   16827   A00336:A00336:HNKCMDRXX:2:2121:18511:18646      7       -       65M6S   GGGCCGGGGACCTCCCTGGTCACACACCTTCTTCCCTAGACACCCCACACTTTGTGTTTCAGACCTA
        bedtools getfasta return does not include suspending seq corresponding to S in cigar

        Note:
        At the end of chr, if the requested number of bases exceeds the boundaries of the chr, bedtools slop will clip the feature accordingly.
        This may cause index error in Convert2XM step.
        e.g. for reads at the end of chr M 
        This step consumes high amount of memory.
        '''
        print "Get genome context."
        
        # reference genome and index used for alignment
        fa = self.fa
        fai = self.fa + '.fai'

        suffix = os.path.basename(self.sam)
        if os.path.exists('bamtobed.{}.script'.format(suffix)):
            print "tempory file for getfasta exist! Quit."
            quit()

        with open('bamtobed.{}.script'.format(suffix), 'w') as output:
            print>>output, '#!/bin/sh'

            # should be the same as using: bedtools slop -s -r 2 -l 0, for both strands
            print>>output, '{} view -h {} -b | bedtools bamtobed -i - -cigar > {}.bed'.format(pathSamtools ,self.sam, suffix)
            command = '{} slop -i {}.bed -g {} -r 2 -l 0 -s | {} getfasta -fi {} -bed - -s -bedOut > {}'.format(pathBedtools, suffix, fai, pathBedtools, fa, 'methylation.temp{}.fa'.format(suffix))
            print>>output, command
        try:    
            check_call('sh bamtobed.{}.script'.format(suffix), shell=True)
        except:
            print "Getfasta using bedtools Error."
            quit()

        # For read from pair-end sam, do not include the readpair information (/1 or /2 for Read1 and Read2) to be consistent with the readID in Sam file. 
        with open('methylation.temp{}.fa'.format(suffix)) as f1:
            key = (line.strip().split('\t')[3].split('/')[0] for line in f1) # remove /1 or /2 information for readID to match the sam file ID.
            with open('methylation.temp{}.fa'.format(suffix)) as f2:
                value = (line.strip().split('\t')[7] for line in f2)
                self.contextgenome = dict(zip(key, value)) # {read ID: context genome ...}
        
        return ('methylation.temp{}.fa'.format(suffix), 'bamtobed.{}.script'.format(suffix), '{}.bed'.format(suffix))
        # the temp files generated: e.g. 'methylation.tempAddtagsam_00.fa', 'bamtobed.Addtagsam_00.script', 'forwardstrand.Addtagsam_00.bed', 'reversestrand.Addtagsam_00.bed'

    def FindMethylation(self, ID, cigar):
        '''
        convert the methylation status in deconvolution report into methylation tag,
        return methylation tag with status for every position in seq, 
        e.g. '......C..C.C..M..CM.......C.CCC..C...CC....',
        the returned methylation tag without consideration of deletion or insertion or suspension in cigar (DIS in cigar);
        the deletion or insertion or suspension will be considered in ConvertSeq step.
        '''
        
        status = self.metylationstatus.get(ID, None) 
        '''
        If deconvolution report has ID and methylation information:
            status = information, e.g. C4C6C8C10C11C22C24C29C35C37C40C41C43C45C46C48C49C51C52C53C55C56C58
        elif deconvolution report has ID but no methylation information:
            status = -1 (defined in MethylationStatus step)
        else: deconvolution report does not have ID
            status = None
        '''
        length = sum([int(re.findall('(\d+)[MIS]', item)[0]) for item in cigar if 'D' not in item])
        # dic_C: saving the methylation status for positions in seq based on deconvolution report
        # For Flag 16, the the numbering of methulation status in report is based on the RCseq of column 10
        tag = [] # tag saving status '...C...M', . -> non C; C -> unmethylation C; M -> methylated C
        if status: # ID exist in deconvolution report
            if status != -1: # ID has the methylation information
                key = [re.findall('[CM](\d+)', item)[0] for item in re.findall('[CM]\d+', status)]
                value = [re.findall('([CM])\d+', item)[0] for item in re.findall('[CM]\d+', status)]
                dic_C = dict(zip(key, value)) # {index: status} e.g. {'51': 'C', '29': 'C, '43': 'M', ..}
                
                for i in range(0, length):
                    if dic_C.has_key(str(i)):
                        tag.append(dic_C[str(i)])
                    else:
                        tag.append('.')
            else: # No C in the original read
                for i in range(0, length):
                    tag.append('.')  
        else: # return [] for ID not included in deconvolution report 
            pass

        return tag # e.g. '......C..C.C..M..CM.......C.CCC..C...CC....', if tag == [] -> ID is not included in deconvolution report 

    def ConvertSeq(self, seq, cigar):
        '''
        Convert col10 seq (or corresponding methylation tag) to sequence with consideration of deletion or insertion or suspension.
        Remove the suspension part at both ends if existing, suspesion S always at the begining or end of the cigar,
        so the suspension part does not have corresponding information in XM tag.

        use '^' to represent bases corresponding to D in cigar.
        use '-' to represent bases corresponding to I in cigar.

        seq = 'GTGATGCAGCTCTTCGGCCTGGTTAACACCCTTCTGGCCAATGACCCAACATCCCTTCGGAAAANCTC'
        cigar = ['60M', '2I', '2M', '1D', '4M']
        return GTGATGCAGCTCTTCGGCCTGGTTAACACCCTTCTGGCCAATGACCCAACATCCCTTCGG--AA^NCTC

        seq = 'ATAGACTACATCTANATAATTGCCTCTTCAAATAGTCGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT'
        cigar = ['36S', '64M']
        methylation = 'C5C8C11C22C23C25C28M36C42C43C44M48C49C50C54C55C56C60C61C62C66C67M68C72C73C74C78C79C80C84M85C86C90C91C92C96C97C98'
        return:
        'CGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT'
        'M.....CCC...MCC...CCC...CCC...CCM...CCC...CCC...CMC...CCC...CCC.'

        seq = 'CCTAACCCTANCCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCNAACCNTAACNCTA'
        cigar = ['41M', '1I', '41M', '1D', '17M']
        methylation = 'C0C1C5C6C7C11C12C13C17C18C19C23C24C25C29C30C31C35C36C37C41C42C43C44C47C48C49C53C54C55C59C60M61C65C66C67C71C72C73C77C78C79C83C84C85C89C90C95C97'
        return:
        'CCTAACCCTANCCCTAACCCTAACCCTAACCCTAACCCTAA-CCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA^CCCNAACCNTAACNCTA'
        'CC...CCC...CCC...CCC...CCC...CCC...CCC...-CCC..CCC...CCC...CCM...CCC...CCC...CCC...^CCC...CC....C.C..'

        For Flag 16, use RC of SEQ and reversed cigar list as input 
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
    
    def FindContext(self, pos, seq):
        '''
        Use to find the methylation context for position based on genome sequence

        pos: index in context genome
        seq: context genome (including the 2bp downstream)
        
        return Z:CpG, H:CHH, X:CHG
        return 'U' means CN or CHN
        return '.' means seq[pos] is not C in genome
        '''
        seq = seq.upper() # genome seq has both upper and lower case, convert to upper, ^ or - will not change
        
        CG = ['CGA', 'CGT', 'CGC','CGG']
        CHH = ['CAA', 'CAT', 'CAC', 'CTT', 'CTA', 'CTC', 'CCC', 'CCA', 'CCT']
        CHG = ['CAG', 'CTG', 'CCG']

        # dic has the same naming as bismark, but use upper letter (do not distinguish methylated or not here)
        dic = dict.fromkeys(CG, 'Z')
        dic.update(dict.fromkeys(CHH, 'H'))
        dic.update(dict.fromkeys(CHG, 'X'))
        
        # return 'U' means having N, while return '.' means seq[pos] is not C
        if seq[pos].upper() == "C":
            return dic.get(seq[pos: pos+3], 'U')
        else: # pos in genome is not C (e.g. -- corresponding to insertion in cigar or mismatch base)
            return '.'

    def Convert2XM(self, ID, tag, seq, cigar):
        '''
        convert the methylation tag based on denvolution report into XM tag.
        
        use tag, seq without addition of deletion or insertion as input:
            tag: methylation tag based on deconvolution report, generated by FindMethylation.
        
        consider the deletion and insertion in this step
        '''
        XM = []

        # include deletion and insertion in the SEQ and context genome and methylation tag
        # newseq = self.ConvertSeq(seq, cigar) # seq with addition of deletion and insertion, e.g. GTGATGCAGCTCTTCGGCCTGGTTAACACCCTTCTGGCCAATGACCCAACATCCCTTCGG--AA^NCTC
        newtag = self.ConvertSeq(self.FindMethylation(ID, cigar), cigar) # methylation tag with addition of deletion and insertion, ....C.C.C.CC..........C.C....C.....C.C..CC.C.CC.CC.CCC.CC.C.--..^....
        # newseq and newtag bases are matching
        
        genome = self.contextgenome[ID] # genome sequence including 2bp downstream
        
        # See BedTools_README for better explaination.
        key = [i for i in range(len(newtag)) if newtag[i] != '-']
        # at the end of chr, bedtools slop will clip the feature accordingly if equested number of bases exceeds the boundaries
        # drop this kind of read if can not get 2bp downstream 
        if len(genome) == len(key)+2:
            # pair the index of newseq, newtag with the index in genome (skip '-' in newseq/newtag)
            dic_index = {} # {index in newtag: index in genome}
            value = [i for i in range(len(genome)-2)]
            dic_index = dict(zip(key, value))
            
            # Find methylation context for each position

            for i in range(len(newtag)): 
                if newtag[i] == '.': # nonC based on deconvolution report, so even the corresponding genome is a C, still add '.' for XM tag.
                    XM.append('.')
                elif newtag[i] == '-':
                    XM.append('.')
                elif newtag[i] == 'C': # unmethylated C based on deconvolution report
                    XM.append(self.FindContext(dic_index[i], genome).lower()) # dic[i]: index in genome
                elif newtag[i] == 'M': # methylated C based on deconvolution report
                    XM.append(self.FindContext(dic_index[i], genome))
                else: # e.g. ^ corresponding to the D in cigar, do not add any symbol
                    pass
            return ''.join(XM)
        else:
            print "Read {} reaches the end for chromsome and can not get methylation context.".format(ID)
            return None


    def CreateXM(self, output_sam, remove_sam):
        '''
        Create XM tag following the bismark principle based on deconvolution report

        Logic:
        Creat XM tag for methylation status corresponing the position in SEQ;
        the suspending part in SEQ (S in cigar) does not have positions in XM tag, so len(XM) = cigar M/I
        
        To report XM tag based on deconvolution report:
            nonC: not a C (or not Methylated C)
            If a position is a C (or Methylated C) in deconvolution report, and also C on genome:
                Find the methylation context in genome, use x,h,z,u or Capital in XM tag.
            If a position is nonC in either deconvolution report or genome:
                Use '.' in XM tag in the position in SEQ
            If a position is a C on genome but is nonC in deconvolution report:
                Use '.' in XM tag.
            If a position is C in deconvolution report, but nonC on genome (mismatch) or this base is a insertion (I in cigar) not present in genome:
                Use '.' in XM tag.

        output_sam: reads (-F 4) having methylation status and added XM tag, having header
        remove_sam : reads not aligned (-f 4) or reads not having methylation status in report, having header
        '''
       
        output = open(output_sam, 'w')
        remove = open(remove_sam, 'w')
        with open(self.sam) as f:
            for line in f:
                if not line.startswith('@'):
                    line = line.strip().split('\t')

                    if not Flag4(line[1]): # -F 4, aligned
                    
                        # check cigar
                        if 'X' in line[5] or 'N' in line[5] or 'H' in line[5] or 'P' in line[5] or '=' in line[5]:
                            print "Cigar containing string not MIDS: {}".format(line[0])
                            continue
                        
                        cigar = re.findall('\d+[MDIS]', line[5]) # '4M1D72M' -> ['4M', '1D', '72M']
                        seq = line[9] # col10 in sam file

                        '''
                        Flag 16:
                        Deconvolution numbering is based on the original fastqseq, which is the RC of SEQ,
                        so use RC of SEQ for methylation status analysis and therefore reverse the cigar too,
                        at the end reverse back, so positions in XM tag corrrepond to SEQ
                        '''
                        if Flag16(line[1]): 
                            seq = RC_seq(seq)
                            cigar = cigar[::-1]
                            tag = self.FindMethylation(line[0], cigar) # tag: methylation tag based on deconvolution report
                            if not tag: # save read for which readID is not included in deconvolution report to 'remove'
                                print>>remove, '\t'.join(line)
                                continue

                            XM = self.Convert2XM(line[0], tag, seq, cigar)
                            if not XM: # at the end of chromosome, can not get 2bp downstream
                                print>>remove, '\t'.join(line)
                                continue

                            line.append('XM:Z:{}'.format(XM[::-1])) # reverse XM tag for Flag 16
                            line.append('XR:Z:CT')
                            line.append('XG:Z:GA')
                        
                        # Flag 0
                        else:
                            tag = self.FindMethylation(line[0], cigar)
                            if not tag:
                                print>>remove, '\t'.join(line)
                                continue
                            
                            XM = self.Convert2XM(line[0], tag, seq, cigar)
                            if not XM: 
                                print>>remove, '\t'.join(line)
                                continue
                            line.append('XM:Z:{}'.format(XM))
                            line.append('XR:Z:CT')
                            line.append('XG:Z:CT')

                        print>>output, '\t'.join(line)
                    
                    else: # -f 4, not aligned reads
                        print>>remove, '\t'.join(line)
                
                else: # header
                    print>>output, line.strip()
                    print>>remove, line.strip()
        output.close()
        remove.close()
        return 1
        
##########
def RunMethylationExtraction(input):
    '''
    input: report, name e.g. Addtagsam_00, 'samfile.header', fa
    run MethylationExtraction for each split file.
    Note:
    The returned sam file Addtagsam_00.output and Addtagsam_00.remove both have header.
    '''
    print "Input sam file: {}".format(input[1])

    # extract readID existing in sam from report, add header to split sam file
    output = open('{}.addheader'.format(input[1]), 'w') # Addtagsam_00.addheader
    with open(input[2]) as f:
        for line in f:
            print>>output, line.strip()

    with open(input[1]) as f:
        key = (line.strip().split('\t')[0] for line in f if not line.startswith('@')) # readID e.g. chr1-99582540
        dic_sam = dict.fromkeys(key)
        f.seek(0)
        for line in f:
            print>>output, line.strip()
    output.close()
    
    # create deconvolution report (name.deconvoluted_5mC) saving the existing ID in sam file
    output = open('{}.deconvoluted_5mC'.format(input[1]), 'w')
    with open(input[0]) as f:
        for line in f:
            if dic_sam.has_key(line.strip().split('\t')[0]):
                print>>output, line.strip()
    output.close()

    # MethylationExtraction(report, sam, fa)
    extract = MethylationExtraction('{}.deconvoluted_5mC'.format(input[1]), '{}.addheader'.format(input[1]), input[3]) # e.g. Addtagsam_00.deconvoluted_5mC
    extract.MethylationStatus()
    tempfile = extract.GetFasta() # prefix for temporary files
    extract.CreateXM('output.{}'.format(input[1]), 'remove.{}'.format(input[1])) 
    return ('output.{}'.format(input[1]), 'remove.{}'.format(input[1])) # e.g. (Addtagsam_00.output, Addtagsam_00.remove)


##-------mainbody
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--name', help='all output files starting with this base file name', dest='prefix', default='addtag', required=False)
    parser.add_argument('--report', help='deconvolution report', dest='report', required=True)
    parser.add_argument('--input', help='input Read1 sam file', dest='input', required=True)
    parser.add_argument('--dir', help='dir path for output file', dest='output_dir', required=False)
    
    parser.add_argument('--thread', help='split sam files and run multiple threads', action='store_true', dest='thread', required=False)
    parser.add_argument('--smp', help='number of thread', dest='smp', type=int, required=False)
    parser.add_argument('--reference', help='reference.fa', dest='fa', required=True)

    parser.add_argument('--path_to_bedtools', help='specify a path to the bedtools executable', dest='pathBedtools', default='bedtools')
    parser.add_argument('--path_to_samtools', help='specify a path to the samtools executable', dest='pathSamtools', default='samtools')
    
    args = parser.parse_args()

    report = GetFilePath(args.report)
    sam = GetFilePath(args.input)
    fa = GetFilePath(args.fa)

    pathBedtools = args.pathBedtools
    pathSamtools = args.pathSamtools

    print "Add XM tag to the mapping of the deconvoluted reads using version: {}".format(__version__)
    print 'Start at {}\n'.format(datetime.datetime.now())

    if args.output_dir:
        dir = args.output_dir
    else: # without --output_dir, saving the output files under working dir
        dir = os.getcwd()

    output_sam = os.path.join(dir, '{}.XMtag.sam'.format(args.prefix)) # .e.g AccTar_deconvolution.XMtag.sam, file saving the reads with methylation information
    remove_sam = os.path.join(dir, '{}.noXMtag.sam'.format(args.prefix)) # .e.g AccTar_deconvolution.noXMtag.sam, file saving the reads without methylation information in deconvolution report

    if args.thread and multipleThreads: # multiple thread by spliting the Read1 sam
        # the number of nodes used for multiple threading n
        if args.smp:
            n = args.smp
        else:
            n = int(pathos.helpers.cpu_count())
        print "Use n={} nodes for pathos multiple threading.".format(n)

        localtime = time.asctime(time.localtime())
        suffix = ''.join(localtime.split()[-2].split(':')) # '151542'
        cwd = os.getcwd()
        os.mkdir('AddtagTmp{}'.format(suffix))
        os.chdir('AddtagTmp{}'.format(suffix))
        print 'Generate a tmp dir {}'.format('AddtagTmp{}'.format(suffix))

        # thread = multiprocessing.cpu_count() # pathos.multiprocessing.cpu_count(), count available thread        
        with open('samtools.script{}'.format(suffix), 'w') as output:
            print>>output, '#!/bin/sh'
            print>>output, '{} view -H {} > samfile.header'.format(pathSamtools, sam)
            print>>output, '{} view {} > samfile.withoutheader.sam'.format(pathSamtools, sam) # split sam file without header, otherwise Addtagsam_00 has duplicated header
        
        try:    
            check_call(['sh', 'samtools.script{}'.format(suffix)], shell=False)
        except:
            print "samtools Error."
            quit()

        # split into subfiles containing 30 million mapping: Addtagsam_00, Addtagsam_01 ...
        # return non-zero exit status 1 'split: output file suffixes exhausted' if needs more than 99 subfiles; so maximum mapping is 2970 million
        try:
            command = 'split -l 30000000 -d -a 2 {} Addtagsam_'.format('samfile.withoutheader.sam')
            check_call(command, shell=True)
            print "Split Read1 sam done."
        except Exception as e: 
            print e
            print 'To solve this error, either splitting the input sam file or changing the number of alignments in each sub file.'
            print 'See README.pdf Addition of XM tag: AddXMtag.py note 6 for details.'
            quit()

        # Use sort to Makesure the split file is in the same order as input: Addtagsam_00, Addtagsam_01 ...
        sam_list_temp = sorted([item for item in os.listdir(os.getcwd()) if re.findall('^Addtagsam_\d\d$', item)])
        sam_list = [(report, item, 'samfile.header', fa) for item in sam_list_temp] # [(path/AccRep1_deconvoluted_5mC.txt, Addtagsam_00, samfile.header) ...]
        # sam_list = [(report, item, 'samfile.header') for item in os.listdir(os.getcwd()) if re.findall('Addtagsam_\d\d$', item)]

        os.remove('samfile.withoutheader.sam')

        print "Create XM tag containing methylation status."
        print "Deconvolution report file: {}\n".format(report)
        pool = ProcessPool(nodes=n) # nodes=4
        ls_output = pool.map(RunMethylationExtraction, sam_list) # ls_output: [(Addtagsam_00.output, Addtagsam_00.remove) ...]
        pool.close()
        pool.join()

        with open(output_sam, 'w') as output:
            with open('samfile.header') as f:
                for line in f:
                    print>>output, line.strip()
            # The returned sam file Addtagsam_00.output and Addtagsam_00.remove both have header, need to remove header for each split
            for item in ls_output:
                with open(item[0]) as f:
                    for line in f:
                        if not line.startswith('@'): # skip header for split output sam
                            print>>output, line.strip()
        
        with open(remove_sam, 'w') as output:
            with open('samfile.header') as f:
                for line in f:
                    print>>output, line.strip()
            for item in ls_output:
                with open(item[1]) as f:
                    for line in f:
                        if not line.startswith('@'):
                            print>>output, line.strip()
        os.chdir(cwd)
        command = ['rm', '-r', 'AddtagTmp{}'.format(suffix)]
        check_call(command, shell=False)

    else: # Single thread, do not split Read1 sam
        print "Create XM tag containing methylation status."
        print "Input sam file: {}".format(sam)
        print "Deconvolution report file: {}".format(report)
        extract = MethylationExtraction(report, sam, fa)
        extract.MethylationStatus()
        tempfile = extract.GetFasta() # temporary files
        extract.CreateXM(output_sam, remove_sam)
        for temp in tempfile:
            os.remove(temp)
    
    print
    print "Output files are:"
    print output_sam
    print remove_sam
    print "Done at {}.".format(datetime.datetime.now())