# -*- coding: utf-8 -*-
# python2.7
'''
# Logs:

Created on Oct 27, 2020
Modifications on Aug 20:2022:
    change check_call to shell=False
Modifications on Feb 23, 2024:
    Modify the CalculateBayesianInference function and CreateQualityscore to address the issue of division by zero;
    I have confirmed this change does not change the counting if sum(Base_[A,T,C,G]) != 0.
'''
'''
Requirement:
(1) Based on python2.7, need pandas and pathos

Purpose:
Based on the Comparison between Conversion_R1 and reference, R2 and reference at each sequencing cycle,
calculate the conditional probability for REF (e.g. E:(R1=G, R2=A), H:(REF=G)) using BayesianInference.
The returning probability in name.BaseCalibration.probability could be used to calibrate the sequencing quality using Deconvolution.py.

Logic and Explanation:
    See class BayesianInference and BayesianDemo.xlsx for details.

Usage:
$python BaseCalibration.py --Read1 TestSeq_deconvoluted_R1.compareBase.txt --Read2 TestSeq_R2.compareBase.txt --name TestSeq

--Read1 and --Read2:
    CompareBase files generated by CompareBase.py using Conversion_R1 and R2.
    e.g.
    ReadID, mapq, cycle, chr, 1-coordination position, Reference_base(REF), Read_base, Read_baseQuality
    A00336:A00336:HNKCMDRXX:1:1104:29044:29684      42      0       chr4    148667385       G       G       37
--name: 
    Output files saving the count and probality matrix and conditional probability at each cycle.
    Output files are saved in dir where the script is run.
    
    name.BaseCalibration.table, e.g.:
    # cycle=0 Read1 count table
            Base_A   Base_C  Base_G  Base_T     sum
    REF_A  714406    24943   17692   24086   781127
    REF_C    7940  1477562    7430   17286  1510218
    REF_G    6468    15062  525000   14206   560736
    REF_T   11975    32346   11446  987051  1042818
    
    name.BaseCalibration.probability, 
    e.g.:
    <Cycle, Read1Base, Read2Base, REF, probability, corresponding quality score>
    Cycle=0 R1=A    R2=T    REF=A   0.674802838049182       &
    Cycle=0 R1=A    R2=T    REF=T   0.3174817528706056      #
    Note:
    Means that under condition R1=A, R2=T, the probability that REF=A is 0.674802838049182 and therefore the converted quality score is '&'.

Note:
(1) storage and memory, current setting smp=nodes=8 for multiple threading
Generate tempory files as large as Read1 and Read2.
Use mulitple threading to speed up (nodes=smp) by splitting based on cycle.
Need enough memory to save data (size around 1/100 of Read1) for pandas at each cycle.
If the input compareBase.txt is 22G, it requires vmem=2.783G, maxvmem=3.320G for running one thread.
(2) Read/Base used for analysis
Currently do NOT apply mapq or sequencing quality cutoff on read or base, means use all the bases in input file for analysis.
Could add filter for mapq or sequencing quality if neccessary.
(3) tempory dir is named based on time.
'''
try:
    import argparse
    import os, sys
    from subprocess import check_call
    from itertools import permutations
    import pandas as pd
    from collections import OrderedDict
    import math
    from pathos.pools import ProcessPool

    script_dir = os.path.dirname(__file__)
    sys.path.append(script_dir)
    from getpath import GetFilePath, CreatePrefix
except:
    print "Import module error."
    quit()


class BayesianInference:
    '''
    Calculate the conditional probability R1, R2 mismatch for REF=A,C,T,G using BayesianInference at certain cycle=n
    The returning probability is used to calibrate the R1, R2 mismatch during deconvolution step.
    '''
    def __init__(self, input_R1, input_R2, Base_R1, Base_R2):
        self.input_R1 = input_R1 # Read1 compareBase at cycle=n
        self.input_R2 = input_R2 # Read1 compareBase at cycle=n
        
        self.dic_type_count_R1 = self.CreateMatrix(self.input_R1)
        self.dic_type_count_R2 = self.CreateMatrix(self.input_R2)
        
        # conditional probability under R1=Base_R1 and R2=Base_R2
        self.df_count_R1, self.df_probability_R1, self.df_count_R2, self.df_probability_R2, self.dic_Posteriorprobability = self.CalculateBayesianInference(Base_R1, Base_R2)

    def CreateMatrix(self, input_file):
        '''
        Calculate the number of counts for each combination of REF-Base in a certain cycle,
        input_file is generated using awk column3==cycle
        return a 4X4 matrix as shown in BayesianDemo.xlsx sheet R1
        '''

        ls_type = list(permutations(['A', 'T', 'C', 'G'], 2)) # 12 combinations
        ls_type.extend([('A','A'), ('T','T'), ('C','C'), ('G','G')]) # (REF, ReadBase) A-A, A-T, A-C, A-G
        
        df = pd.read_csv(input_file, header=None, sep='\t', \
        names=['ReadID', 'mapq', 'cycle', 'chr', 'position', 'REF', 'RB', 'RBQuality']) # no header
        
        dic_type_count = dict.fromkeys(ls_type, 0) # {(REF, Base): count}
        for key in ls_type: # total 16 items
            df_sub = df.loc[(df['REF']==key[0]) & (df['RB']==key[1]), :] # e.g. REF=A, Base=A
            dic_type_count[key] = df_sub.shape[0]
        return dic_type_count # sorted(dic_type_count.keys())

    def CreateProblityMatrix(self, dic_type_count):
        '''
        Convert the count matrix in the dic_type_count dictionary to probability matrix
        '''
        data = OrderedDict()
        data['REF_A'] = [dic_type_count[('A', Base)] for Base in ['A', 'C', 'G', 'T']]
        data['REF_C'] = [dic_type_count[('C', Base)] for Base in ['A', 'C', 'G', 'T']]
        data['REF_G'] = [dic_type_count[('G', Base)] for Base in ['A', 'C', 'G', 'T']]
        data['REF_T'] = [dic_type_count[('T', Base)] for Base in ['A', 'C', 'G', 'T']]
        df = pd.DataFrame.from_dict(data, orient='index', columns=['Base_A', 'Base_C', 'Base_G', 'Base_T'])

        df['sum'] = df.sum(axis=1)
        # calculate the probability matrix, e.g. probability for REF=A,Base=A: count_R1_REFA_BaseA / count_R1_REFA_total
        df_probability = df.loc[:,"A":"T"].div(df['sum'], axis=0)

        return df, df_probability

    def CalculateBayesianInference(self, Base_R1, Base_R2):
        '''
        E:(R1=C,R2=A,T,G), H:(REF=A,T,C,G)
        If R1=G and R2!=A (Base_R1=G, Base_R2=A,T,G), calculate the probability that REF=A,C,T,G;

        return a dictionary dic_Posteriorprobability saving all the probability.

        See BayesianDemo.xlsx for details

        # training data used in BayesianDemo.xlsx    
        dic_type_count_R1 = {('C', 'G'):4181, ('A', 'T'): 10592, ('T', 'T'): 4923851, ('T', 'A'): 9121, ('C', 'A'): 5278, ('C', 'T'): 6878, \
        ('G', 'G'): 2675865, ('T', 'C'): 43350, ('G', 'A'): 4792, ('G', 'T'): 8593, \
        ('A', 'G'): 44562, ('G', 'C'): 3494, ('C', 'C'): 3076616, ('T', 'G'): 5565, ('A', 'A'): 4897824, ('A', 'C'): 4707}
        dic_type_count_R2 = {('C', 'G'): 4430, ('A', 'T'): 14671, ('T', 'T'): 5893820, ('T', 'A'): 16697, ('C', 'A'): 9642, ('C', 'T'): 6728, \
        ('G', 'G'): 3166396, ('T', 'C'): 52354, ('G', 'A'): 105364, ('G', 'T'): 8075, \
        ('A', 'G'): 6770, ('G', 'C'): 5828, ('C', 'C'): 3763024, ('T', 'G'): 6269, ('A', 'A'): 5949181, ('A', 'C'): 7428}

        Note:
        at cycle=99, sometimes there is no Base=G e.g.
                Base_A  Base_C  Base_G  Base_T    sum
        REF_A   75114    3695       0     119  78928
        REF_C    1956   40876       0      86  42918
        REF_G    2398    2009       0      62   4469
        REF_T    3354    4208       0   71011  78573
        This will raise error when calculate the BayesianInference. In this case, return an empty {} as dic_Posteriorprobability.
        '''
        df_count_R1, df_probability_R1 = self.CreateProblityMatrix(self.dic_type_count_R1)
        df_count_R2, df_probability_R2 = self.CreateProblityMatrix(self.dic_type_count_R2) # a certain cell value: df_probability_R1.loc['REF_A','Base_G']

        try:
            # Likelihood: P(E|H=A) = P(REF=A, R1=G)*P(REF=A, R2=A), P(E|H=C) = P(REF=C, R1=G)*P(REF=C, R2=A)
            # e.g. R1=G, R2=A conditional probability: df_probability_R1.loc[REF, 'Base_G']*df_probability_R2.loc[REF, 'Base_A']
            # {'A': P(E|H=A), 'C': P(E|H=C), 'T': P(E|H=T), 'G': P(E|H=G)}
            dic_Likelihood = dict(list(zip(['A', 'C', 'G', 'T'], \
            [df_probability_R1.loc[REF, 'Base_{}'.format(Base_R1)]*df_probability_R2.loc[REF, 'Base_{}'.format(Base_R2)] for REF in ['REF_A', 'REF_C', 'REF_G', 'REF_T']])))

            # Prior probability: P(H=A) = (count_R1_REF=A + count_R2_REF=A)/(count_R1_REF=All + count_R2_REF=All)
            dic_Priorprobability = dict(list(zip(['A', 'C', 'G', 'T'], \
            [float(df_count_R1.loc[REF, 'sum']+df_count_R2.loc[REF, 'sum'])/(sum(df_count_R1.loc[:,'sum'])+sum(df_count_R2.loc[:,'sum'])) for REF in ['REF_A', 'REF_C', 'REF_G', 'REF_T']])))

            # Evidence: P(E) = P(E|H=A)*P(H=A) + P(E|H=C)*P(H=c) + P(E|H=T)*P(H=T) + P(E|H=G)*P(H=G)
            Evidence = sum([dic_Likelihood[key]*dic_Priorprobability[key] for key in ['A', 'C', 'G', 'T']])
            if Evidence == 0:
                return df_count_R1, df_probability_R1, df_count_R2, df_probability_R2, {}
            else:
                # Posterior Probability: P(H=A|E) = P(E|H=A)*P(H=A)/P(E)
                dic_Posteriorprobability = dict(list(zip(['A', 'C', 'G', 'T'], \
                [dic_Likelihood[key]*dic_Priorprobability[key]/Evidence for key in ['A', 'C', 'G', 'T']])))
                
                return df_count_R1, df_probability_R1, df_count_R2, df_probability_R2, dic_Posteriorprobability
                # dic_Posteriorprobability: {'A': 0.402919587687, 'C': 0.00156303537087, 'T': 0.00448675310282, 'G': 0.591030623839}, means probability that REF=A is 0.402919587687 
        except: # e.g. in case sum(df_count_R1.loc[:,'sum'])+sum(df_count_R2.loc[:,'sum'] == 0
            return df_count_R1, df_probability_R1, df_count_R2, df_probability_R2, {}
        
def CreateQualityscore(probability):
    '''
    probability=0.6252511082206418, return &, means QualityScore=5
    probability=0.34212911890837955, return #, means QualityScore=2
    probability=1, return 'NA'
    '''
    try:
        return chr(int(math.ceil(-10*math.log10(1-probability)))+33)
    except: # in case probability==1
        return 'NA'

def Calculateprobability(i):
    '''
    i is the cycle
    Calculate all the probability for all the combinations=48 in each cycle.
    output files are probability and table for cycle i
    '''
    output_table = open('table.cycle{}'.format(i), 'w')
    output_probability = open('probability.cycle{}'.format(i), 'w')

    ls_type = list(permutations(['A', 'T', 'C', 'G'], 2)) # 12 combinations for R1 R2 mismatch
    for item in ls_type:
        # 'R1.cycle{}'.format(i) and 'R2.cycle{}'.format(i) is the split file for cycle i
        temp = BayesianInference('R1.cycle{}'.format(i), 'R2.cycle{}'.format(i), item[0], item[1]) # R1=item[0], R2=item[1]
        
        # print df_count_R1 and df_count_R2 at each cycle, only need to print once for each cycle since the count should be the same
        if item == ('G', 'A'): 
            print>>output_table, '# cycle={} Read1 count table'.format(i)
            print>>output_table, temp.df_count_R1
            print>>output_table, '# cycle={} Read2 count table'.format(i)
            print>>output_table, temp.df_count_R2

        # print probability and save it in returning dic_probability
        for key in ['A', 'T', 'C', 'G']:
            ls_key = ['Cycle={}'.format(i), 'R1={}'.format(item[0]), 'R2={}'.format(item[1]), 'REF={}'.format(key)]
            if temp.dic_Posteriorprobability: # key: ['A', 'C', 'T', 'G']
                ls_value = [str(temp.dic_Posteriorprobability[key]), CreateQualityscore(temp.dic_Posteriorprobability[key])]
            else: # probability not available, do not save in the dic_probability
                ls_value = ['NA', 'NA']
            ls_key.extend(ls_value)

            print>>output_probability, '\t'.join(ls_key)
            # Cycle=1   R1=G    R2=A    REF=G   0.591030623839  %
            # when R1=G, R2=A, the probability of REF=G is 0.591030623839, correponding to quality score %
    output_table.close()
    output_probability.close()
    return 'table.cycle{}'.format(i), 'probability.cycle{}'.format(i) # the name of output files for this cycle
    

def main(Read1, Read2, name):
    '''
    Split the large Read1 and Read2 file based on the cycle
    Run Calculateprobability() for each cycle using multiple threads
    '''
    prefix = CreatePrefix()
    
    Read1 = GetFilePath(Read1)
    Read2 = GetFilePath(Read2)

    CWD = os.getcwd()

    os.mkdir('BaseCalibration{}'.format(prefix)) # BaseCalibration172545
    os.chdir('BaseCalibration{}'.format(prefix))

    # split based on the cycle, note this step will generate tempory files as large as Read1 and Read2
    command = 'awk -v OFS=\'\\t\' \'{{print >"{}"$3}}\' {}'.format('R1.cycle', Read1) # e.g. R1.cycle0, R1.cycle1 ... R1.cycle99
    check_call(command, shell=True)
    command = 'awk -v OFS=\'\\t\' \'{{print >"{}"$3}}\' {}'.format('R2.cycle', Read2) # e.g. R2.cycle0, R2.cycle1 ... R2.cycle99
    check_call(command, shell=True)

    # maximum illumina sequencing cycle in the fastq files
    ls_cycle = []
    for item in os.listdir(os.getcwd()):
        if 'R1.cycle' in item:
            ls_cycle.append(int(item.replace('R1.cycle', '')))
    cycle = max(ls_cycle)

    # Calculate all the probability for all the combinations=48 in each cycle.
    # Using multiple thread
    pool = ProcessPool() # nodes=ncpu
    ls_return = pool.map(Calculateprobability, range(0, cycle+1))
    # [('table.cycle1', 'probability.cycle1'), ('table.cycle2', 'probability.cycle2'), ...]
    pool.close()
    pool.join()

    output_probability = open(os.path.join(CWD, '{}.BaseCalibration.probability'.format(name)), 'w')
    output_table = open(os.path.join(CWD, '{}.BaseCalibration.table'.format(name)), 'w')

    for item in ls_return:
        with open(item[0]) as f:
            for line in f:
                print>>output_table, line.strip()
        with open(item[1]) as f:
            for line in f:
                print>>output_probability, line.strip()

    output_probability.close()
    output_table.close()

    print "The output files are:"
    print os.path.join(CWD, '{}.BaseCalibration.probability'.format(name))
    print os.path.join(CWD, '{}.BaseCalibration.table'.format(name))

    os.chdir(CWD)
    command = ['rm', '-r', 'BaseCalibration{}'.format(prefix)]
    check_call(command, shell=False)

    return os.path.join(CWD, '{}.BaseCalibration.probability'.format(name)), os.path.join(CWD, '{}.BaseCalibration.table'.format(name))
    
##-------mainbody
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Read1', help='Deconvolution_R1 compareBase file', dest='Read1')
    parser.add_argument('--Read2', help='Read compareBase file', dest='Read2')
    parser.add_argument('--name', help='name for output file', dest='name')
    args = parser.parse_args()

    main(args.Read1, args.Read2, args.name)


