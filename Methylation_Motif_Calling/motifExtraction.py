# python3.7
'''
# Logs:

Created on Feb 2, 2022
'''
'''
Used to extract a certain size of sequence (e.g. 8bp) containing a methylated C or C at a given position based on Methyl-SNP-seq reads.
The extracted sequences can be used for methylation motif calling using identifyMotif.py and clusterMotif.py.

Logic:
If --contigs is provided:
    Use the reads mapped to the target contig(s) for sequenc extraction
Based on the methylation report, extract the nearby sequence of methylated C (methylC.sequence) or unmethylated C (unmethylC.sequence) or all C (allC.sequence).

Usage:
$python motifExtraction.py --input Ecoli_ABK_Deconvolution_R1.fq --report Ecoli_ABK.Deconvolution.5mC\
    --sam Ecoli_ABK.sam --name Ecoli_ABK_Node-27-70 \
    --contigs NODE_27_length_74030_cov_17.136755 NODE_70_length_61531_cov_18.801319 -left 1 -right 4

--input: Methyl-SNP-seq Deconvolution Read fastq file
--report: corresponding Deconvolution methylation report
--name: name for all output files, Default TestMotif
    Generate Output files:
    name.Deconvolution_R1.fq: containing reads used for seq extraction if --contigs is on
    name.Deconvolution.5mC: containing methylation information for reads in name.Deconvolution_R1.fq
    name_motif_methylC.sequence, name_motif_unmethylC.sequence and name_motif_allC.sequence, e.g.
    e.g.
    CCAGGC
    CCAGGA
    GCAACA
    ..
    Note:
    name_motif_methylC.sequence and name_motif_allC.sequence can be used as input files for identifyMotif.py.

--left/-l INT default 1, --right/-r INT default 4:
    -l and -r define the number of bases upstream of downstream of a methylated or unmethylated C. 
    With default setting, seq is x5mCxxxx in methylC.sequence; xCxxxx in unmethylC.sequence; xC/5mCxxxx in allC.sequence.

--contigs: chr(s), nargs="+",  Not required
    If provided, only the reads mapped to the given chr(s) are used to extract sequences.
    These given chr(s) must be present in the sam file otherwise no read will be used.
--sam: Not required
    Deconvolution Reads mapped to genome assembly, only required if --contigs is provided.

'''
try:
    import re
    import regex
    import argparse
except:
    print("Module Error")

class subread:
    '''
    Use to extract the reads and corresponding methylation information mapped to the given contigs/chrs
    input and report: Deconvolution_R1.fq and corresponding Deconvolution.5mC methylation report
    contigs: chr existing in sam file based on which reads are clustered

    Output files:
    prefix.Deconvolution_R1.fq
    prefix.Deconvolution.5mC
    '''
    def __init__(self, input, report, sam, contigs, prefix):
        self.report = report # deconvolution report
        self.input = input # Deconvolution.fq
        self.sam = sam # sam file, Deconvolution Reads mapped to genome assembly
        self.prefix = prefix
        self.dic_ABK = dict.fromkeys(contigs) # contigs is a list containing the chr(s), ['NODE_27_length_74030_cov_17.136755', 'NODE_70_length_61531_cov_18.801319']

        dic_read = self.ExtractReadID() # a dic having the ReadID in sam mapped to target contigs.
        if dic_read:
            self.ExtractReport(dic_read)
            self.ExtractRead(dic_read)
        else:
            print("No reads mapped to the given contigs. Quit.")
            quit()

    def ExtractReadID(self):
        '''
        Extract the entries in sam file having chr=contigs
        '''
        ls_read = []
        with open(self.sam) as f:
            for line in f:
                if not line.startswith('@'):
                    if self.dic_ABK.has_key(line.strip().split()[2]):
                        ls_read.append(line.strip().split()[0])
        return dict.fromkeys(ls_read)

    def ExtractReport(self, dic_read):
        '''
        Extract the methylation report for Deconvolution Read existing in the dic
        '''
        output = open('{}.Deconvolution.5mC'.format(self.prefix), 'w')
        with open(self.report) as f:
            for line in f:
                if dic_read.has_key(line.strip().split()[0]):
                    print(line.strip(), end='\n', file=output)
        output.close()

    def ExtractRead(self, dic_read):
        '''
        Extract Deconvolution Read existing in the dic
        '''
        output = open('{}.Deconvolution_R1.fq'.format(self.prefix), 'w')
        with open(self.input) as f:
            line1, line2, line3, line4 = f.readline().strip(), f.readline().strip(), f.readline().strip(), f.readline().strip()
            while line1:
                if dic_read.has_key(line1.strip().split()[0].replace('@', '')):
                    print(line1, end='\n', file=output)
                    print(line2, end='\n', file=output)
                    print(line3, end='\n', file=output)
                    print(line4, end='\n', file=output)
                line1, line2, line3, line4 = f.readline().strip(), f.readline().strip(), f.readline().strip(), f.readline().strip()
        output.close()



class motif:
    '''
    extract the motif: x5mCxxxx -> l=1, r=4, total=6 bp motif including the 5mC at +1 position
    Output files: prefix_motif_methylC.sequence and prefix_motif_unmethylC.sequence
    '''
    def __init__(self, input, report, prefix, l, r):
        self.report = report # deconvolution report
        self.input = input
        self.prefix = prefix
        self.l, self.r = l, r

        self.ExtractMotifForMethylC()
        self.ExtractMotifForUnMethylC()
        self.ExtractMotifForAllC()

    def ExtractMotifForMethylC(self):
        '''
        report: Deconvolution.5mC
        input: Deconvolution_Calibration_R1.fq
        x5mCxxxx -> l=1, r=4, total=6 bp motif including the 5mC at +1 position
        '''
        count = 0
        output_file = '{}_motif_methylC.sequence'.format(self.prefix)
        output = open(output_file, 'w')
        with open(self.report) as f1:
            with open(self.input) as f2:
                line1 = f2.readline().strip()
                line2 = f2.readline().strip()
                line3 = f2.readline().strip()
                line4 = f2.readline().strip()
                line_5mC = f1.readline().strip()
                while line1:
                    # A00336:A00336:HV7F7DRXX:1:1101:10004:10019, @A00336:A00336:HV7F7DRXX:1:1101:10004:10019 1:N:0:ATTACTCGCCTATCCT
                    try: # if No C in Deconvolution.5mC report
                        if line1.split()[0].replace('@', '') == line_5mC.split()[0]: # same ID
                            if next(regex.finditer('M(\d+)', line_5mC.split()[1], overlapped=False), None): # read has methylated C
                                methyl_index = [int(index) for index in re.findall('M(\d+)', line_5mC.split()[1])]
                                for index in methyl_index:
                                    if index-self.l>=0 and index+self.r<=len(line2)-1: # -5bp C +5bp, 11bp motif
                                        print(line2[index-self.l:index+self.r+1], end='\n', file=output)
                                        count += 1
                    except:
                        pass
                    line1 = f2.readline().strip()
                    line2 = f2.readline().strip()
                    line3 = f2.readline().strip()
                    line4 = f2.readline().strip()
                    line_5mC = f1.readline().strip()
        output.close()
        print("There are {} kmers having methylated C.".format(count))

    def ExtractMotifForUnMethylC(self):
        '''
        file1: Deconvolution.5mC
        file2: Deconvolution_Calibration_R1.fq.gz
        x5mCxxxx -> l=1, r=4, total=6 bp motif including the 5mC at +1 position
        '''
        count = 0
        output_file = '{}_motif_unmethylC.sequence'.format(self.prefix)
        output = open(output_file, 'w')
        with open(self.report) as f1:
            with open(self.input) as f2:
                line1 = f2.readline().strip()
                line2 = f2.readline().strip()
                line3 = f2.readline().strip()
                line4 = f2.readline().strip()
                line_5mC = f1.readline().strip()
                while line1:
                    # A00336:A00336:HV7F7DRXX:1:1101:10004:10019, @A00336:A00336:HV7F7DRXX:1:1101:10004:10019 1:N:0:ATTACTCGCCTATCCT
                    try: # if No C in Deconvolution.5mC report
                        if line1.split()[0].replace('@', '') == line_5mC.split()[0]: # same ID
                            if next(regex.finditer('C(\d+)', line_5mC.split()[1], overlapped=False), None): # read has unmethylated C
                                methyl_index = [int(index) for index in re.findall('C(\d+)', line_5mC.split()[1])]
                                for index in methyl_index:
                                    if index-self.l>=0 and index+self.r<=len(line2)-1: # -5bp C +5bp, 11bp motif
                                        print(line2[index-self.l:index+self.r+1], end='\n', file=output)
                                        count += 1
                    except:
                        pass
                    line1 = f2.readline().strip()
                    line2 = f2.readline().strip()
                    line3 = f2.readline().strip()
                    line4 = f2.readline().strip()
                    line_5mC = f1.readline().strip()
        output.close()
        print("There are {} kmers having unmethylated C.".format(count))

    def ExtractMotifForAllC(self):
        '''
        file1: Deconvolution.5mC
        file2: Deconvolution_Calibration_R1.fq.gz
        x5mCxxxx -> l=1, r=4, total=6 bp motif including the 5mC at +1 position
        '''
        count = 0
        output_file = '{}_motif_allC.sequence'.format(self.prefix)
        output = open(output_file, 'w')
        with open(self.report) as f1:
            with open(self.input) as f2:
                line1 = f2.readline().strip()
                line2 = f2.readline().strip()
                line3 = f2.readline().strip()
                line4 = f2.readline().strip()
                line_5mC = f1.readline().strip()
                while line1:
                    # A00336:A00336:HV7F7DRXX:1:1101:10004:10019, @A00336:A00336:HV7F7DRXX:1:1101:10004:10019 1:N:0:ATTACTCGCCTATCCT
                    try: # if No C in Deconvolution.5mC report
                        if line1.split()[0].replace('@', '') == line_5mC.split()[0]: # same ID
                            if next(regex.finditer('[M|C](\d+)', line_5mC.split()[1], overlapped=False), None): # read has C
                                methyl_index = [int(index) for index in re.findall('[M|C](\d+)', line_5mC.split()[1])]
                                for index in methyl_index:
                                    if index-self.l>=0 and index+self.r<=len(line2)-1: # -5bp C +5bp, 11bp motif
                                        print(line2[index-self.l:index+self.r+1], end='\n', file=output)
                                        count += 1
                    except:
                        pass
                    line1 = f2.readline().strip()
                    line2 = f2.readline().strip()
                    line3 = f2.readline().strip()
                    line4 = f2.readline().strip()
                    line_5mC = f1.readline().strip()
        output.close()
        print("There are {} kmers having C.".format(count))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', help='deconvolution.5mC report', dest='report')
    parser.add_argument('--input', help='deconvolution Read fq', dest='input')
    parser.add_argument('--name', help='all output files starting with this base file name', dest='prefix', default='TestMotif')

    parser.add_argument('--sam', help='sam file mapping to genome assembly', dest='sam', required=False)
    parser.add_argument('--contigs', nargs="+", help='contigs/chrs used to cluster the reads', dest='contigs', required=False)

    parser.add_argument('-l', '--left', help='The number of base pairs before the 5mC', dest='l', type=int, default=1)
    parser.add_argument('-r', '--right', help='The number of base pairs after the 5mC', dest='r', type=int, default=4) # x5mCxxxx


    args = parser.parse_args()

    # generate: prefix.Deconvolution_R1.fq, prefix.Deconvolution.5mC, which have the reads mapped to contigs
    if args.contigs:
        if args.sam:
            subread(args.input, args.report, args.sam, args.contigs, args.prefix)
            input_fq = '{}.Deconvolution_R1.fq'.format(args.prefix)
            input_report = '{}.Deconvolution.5mC'.format(args.prefix)
        else:
            print('Need to provide a sam file containing the mapping of Deconvoluted reads to genome assembly.')
            quit()

    else: # if no contigs provided, use all the reads in args.input to extract motif region
        input_fq = args.input
        input_report = args.report

    # generate: prefix_motif_methylC.sequence and prefix_motif_unmethylC.sequence
    motif(input_fq, input_report, args.prefix, args.l, args.r)