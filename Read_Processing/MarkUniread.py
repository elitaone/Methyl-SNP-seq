#!/mnt/home/ettwiller/yan/exe/miniconda2/envs/my-python2-env/bin/python
'''
Log:

Created on Aug 16, 2020
Modifications on Aug 26, 2020:
    Add option to filter based on MAPQ
'''
'''
python 2.7, need numpy installed

Use to extract the uniread from input sam file.

Note:
Designed for single-end bowtie2/bwa mem mapping, for which output sam contains Flag 4 labeling unmapped reads.
Do not change the order of entries in the output.

Logic:
    Defined as reliable (uniquely) mapping and save in output if:
    (1) Flag !=4, !=256, !=2048
    (2) XS not in entry or XS in and AS!=XS

Defination of uniquely mapping by aligner:
    bowtie2 default mapping uniread: mapped reads (Flag!=4 and Flag!=256 and Flag!=2048) with only AS tag or AS!=XS tag.
    bwa mem default mapping uniread: mapped reads (Flag!=4 and Flag!=256 and Flag!=2048) with AS!=XS tag.
    Chose to apply MAPQ score or not.

Usage:
$python MarkUniread.py --input sam --output sam --MAPQ default 0

--input:
    a sorted by cooridnation sam
--output:
    a sorted by cooridnation sam
'''
try:
    import re
    import argparse
    import numpy as np
except:
    print "module error"
    quit()
    

def Flag256(Flag):
    '''
    return True if Flag contains Flag 256, Flase if not
    '''
    if len(str(np.binary_repr(int(Flag))))<9: # Flag < 64
        return False
    else:
        return str(np.binary_repr(int(Flag)))[-9] == '1'

def Flag2048(Flag):
    '''
    return True if Flag contains Flag 2048, Flase if not
    '''
    if len(str(np.binary_repr(int(Flag))))<12: # Flag < 64
        return False
    else:
        return str(np.binary_repr(int(Flag)))[-12] == '1'

def Uniread(input_sam, output_sam, mapq):
    '''
    Use to save the uniread in output sam file.
    bowtie2 default mapping uniread: mapped reads (Flag!=4 and Flag!=256) with only AS tag or AS!=XS tag.
    bwa mem default mapping uniread: mapped reads (Flag!=4 and Flag!=256) with AS!=XS tag.
    '''
    total_output = 0
    output = open(output_sam, 'w')
    with open(input_sam) as f:
        line = f.readline().strip()
        while line.startswith('@'): # header
            print>>output, line
            line = f.readline().strip()
        
        if line.split()[1] != '4' and not Flag256(line.split()[1]) and not Flag2048(line.split()[1]) and int(line.split()[4]) >= mapq:
            if 'XS' not in line:
                print>>output, line
                total_output +=1
            else:
                if re.findall('AS:i:(.*?)\s', line)[0] != re.findall('XS:i:(.*?)\s', line)[0]:
                    print>>output, line
                    total_output +=1
        
        for line in f:
            line = line.strip()
            if line.split()[1] != '4' and not Flag256(line.split()[1]) and not Flag2048(line.split()[1]) and int(line.split()[4]) >= mapq:
                if 'XS' not in line:
                    print>>output, line
                    total_output +=1
                else:
                    if re.findall('AS:i:(.*?)\s', line)[0] != re.findall('XS:i:(.*?)\s', line)[0]:
                        print>>output, line
                        total_output +=1

    output.close()
    print "Extract uniquely mapping: with only AS tag or AS!=XS tag."
    print "Input file: {}".format(input_sam)
    print "There are {} number of reads in the output.".format(total_output)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input sam', dest='input')
    parser.add_argument('--output', help='output sam', dest='output')
    parser.add_argument('--MAPQ', help='mapq cutoff', type = int, dest='mapq', default = 0)

    args = parser.parse_args()

    Uniread(args.input, args.output, args.mapq)
