# python3.7, scipy>=1.6.0
'''
Log:
Created on March 3, 2022
Modified on May 4, 2022:
    Add --mode options to control the calculation of Pi_0
'''
'''
Use to identify the significantly methylated motif by comparing the counts in sample and in reference based on Binomial test and Bonferroni correction.
null hypothesis: the methylation level of each motif is not significantly higher than Pi_0 (defined by --mode).
See Ecoli_ABK_test_NODE_38_motif.xlsx for fomula and examples.

Usage:
$python identifyMotif.py --sample Ecoli_motif_methylC.sequence --reference Ecoli_motif_allC.sequence \
    --output Ecoli_SignificantMotif.txt

--sample, --reference:
motif sequences served as positive (seq containing methylated C) and background (seq containing all C) for binomial test.
Generated by motifExtraction.py.
e.g.
CCAGGC
CCAGGA
GCAACA
...

--output: the first line is header
Name    Sample  Reference       pvalue  corrected_alpha Significance
GCCCAGGT        157     285     0.0     6.105378838756945e-09   True
CACCAGGC        282     429     0.0     6.105378838756945e-09   True

--alpha: Default 0.0001
alpha value used to determine the corrected pvalue for significance. 

--cutoff: Default 2
To avoid calling significance of motif seq showing accidentally in sample, 
only the motif seq having count>=cutoff in the reference will be called as significant if adjusted pvalue below sigificance value.

--mode: options [average, top, a float number]
controling the Pi_0 used for binomial test
null hypothesis: the methylation level of seq in sample is not significantly higher than Pi_0
    average: Default mode, Pi_0 = (number of methylated seq in sample) / (number of all seq in reference)
    top: Pi_0 is the 95 percentile of the methylation level of all the seq in sample
    a float value: giving a custom value betweeen 0 and 1 for Pi_0, e.g. 0.67
'''
import pandas as pd
from scipy.stats import binom # scipy 1.8.0
import numpy as np
import argparse

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

class TestSignificance:
    '''
    Find the significant enriched motif seq
    '''
    def __init__(self, sample, reference, alpha, mode):
        self.input_sample = sample
        self.input_reference = reference
        self.alpha = alpha
        self.mode = mode
        
    def CountSeq(self, input_file):
        '''
        count the number of motif in positive and background file
        '''
        dic = {} # {motif: count}
        with open(input_file) as f:
            for line in f:
                if line.strip() in dic:
                    dic[line.strip()] += 1
                else:
                    dic[line.strip()] = 1
        return dic
    
    def CreateDateframe(self):
        '''
        create a dataframe: Name, Count in Sample, Count in Reference
        '''
        dic_sample = self.CountSeq(self.input_sample)
        dic_reference = self.CountSeq(self.input_reference)

        ls_key = [key for key in dict(dic_sample,**dic_reference)]
        ls_sample = [dic_sample.get(key, 0) for key in ls_key]
        ls_reference = [dic_reference.get(key, 0) for key in ls_key]
        ls_methylation = [float(ls_sample[i])/ls_reference[i] if ls_reference[i] > 0 else 0 for i in range(0, len(ls_sample))]

        df = pd.DataFrame({'Name': pd.Series(ls_key), 'Sample': pd.Series(ls_sample), \
            'Reference': pd.Series(ls_reference), 'Methylation': pd.Series(ls_methylation)})
        return df
    
    def BinomTest(self):
        df = self.CreateDateframe()
        
        sum_sample, sum_reference = sum(df['Sample']), sum(df['Reference'])
        if self.mode == 'average': # Pi_0 is the average methylation level
            Pi_0 = sum_sample/sum_reference
        elif self.mode == 'top': # Pi_0 is the 95th of the methylation level
            Pi_0 = np.quantile(df['Methylation'], 0.95)
        elif isfloat(self.mode): # a custom Pi_0, e.g. 0.67
            Pi_0 = float(self.mode)
        print('Pi_0 used for binomial test: {}'.format(Pi_0))
        
        df['pvalue'] = 1 - binom.cdf(df['Sample'], df['Reference'], Pi_0)
        # binom.cdf(k,n,p) return the chances of getting k successes or less in n trials where the success probability is p
        df['corrected_alpha'] = self.alpha/len(df['Sample'])
        df['Significance'] = df['pvalue']<=df['corrected_alpha']
        return df

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--sample', help='positive motif seq', dest='sample')
    parser.add_argument('--reference', help='background motif seq', dest='reference')
    parser.add_argument('--output', help='output file', dest='output')
    parser.add_argument('--alpha', help='alpha value to determine the significance', dest='alpha', type=float, default=0.0001)
    parser.add_argument('--cutoff', help='cutoff of number of motif seq in the reference', dest='cutoff', type=int, default=2)
    parser.add_argument('--mode', help='options for Pi_0 calculation', dest='mode', default='average')

    args = parser.parse_args()

    if args.mode == 'average' or args.mode == 'top':
        pass
    elif isfloat(args.mode):
        if float(args.mode)>=0 and float(args.mode)<=1:
            pass
        else:
            print('The custom Pi_0 needs to be between 0 and 1.')
    else:
        print('Mode for calculation of Pi_0 used for binomial test has to be one of: average, top or a number between 0 and 1.')
        quit()
        
    Result = TestSignificance(args.sample, args.reference, args.alpha, args.mode).BinomTest()

    motif = Result.loc[(Result['Significance']==True) & (Result['Reference']>=args.cutoff)]
    motif.to_csv(args.output, index=False, sep='\t', header=True)

