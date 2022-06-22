# python2.7
'''
Based on python 2.7

Use to remove PCR duplicates for Single end mapping (Deconvolution Read mapping).
Duplicates are defined as reads mapped to the same chr and locus and have the same sequence in col 10.
For duplicates, save the one with highest MAPQ score.

Note:
Here do not consider AS and XS tag, therefore should apply MarkUniread.py before this step if filtering multiple mapping is required.

Usage:
$python MarkUniread.py --input sam --output sam

--input:
    sam file sorted by cooridnation
--output:
    sam file sorted by cooridnation with PCR duplicates removed
'''

from collections import deque
import operator
import argparse

def RemoveDup(ls):
    '''
    Duplicates are defined as reads mapped to the same chr and locus and have the same sequence in col 10.
    Return a list containing the ones with the highest MAPQ if there are duplicates or the only one if there is no duplicate
    '''
    dic = {} # {seq: [index in ls, ...]}, value is a list containing the indexes in ls corresponding to the duplicates
    for i in range(0, len(ls)):
        if dic.has_key(ls[i].split()[9]):
            dic[ls[i].split()[9]].append(i)
        else:
            dic[ls[i].split()[9]] = [i]
    
    dic_return = {} # {seq: index for return entry with higest MAPQ},
    for item in dic:
        if len(dic[item]) > 1:
            dic_mapq = dict(zip(dic[item], [int(ls[i].split()[4]) for i in dic[item]])) # for dplicates: {index in ls: mapq}
            dic_return[item] = max(dic_mapq.iteritems(), key=operator.itemgetter(1))[0] # get the maix MAPQ
        else: # if no duplicate, return the only one
            dic_return[item] = dic[item][0]
            
    return [ls[dic_return[item]] for item in dic_return]

def MarkDup(input_sam, output_sam):
    '''
    Remove the duplicates from input sam file.
    Save the highest MAPQ for duplicates in the output sam.
    '''
    output = open(output_sam, 'w')
    with open(input_sam) as f:
            line = f.readline().strip()
            while line.startswith('@'): # header
                print>>output, line
                line = f.readline().strip()

            q = deque(maxlen=2)
            q.append(line) # first line
            ls = [line]
            line = f.readline().strip() # second line
            while line:
                q.append(line)
                if q[0].split()[2] == q[1].split()[2] and q[0].split()[3] == q[1].split()[3]: # chr and position the same
                    ls.append(q[1])
                else: 
                    print>>output, '\n'.join(RemoveDup(ls))
                    ls = [line]
                line = f.readline().strip() # second line
            
            print>>output, '\n'.join(RemoveDup(ls))

    print "Remove PCR duplicates from input: {}.".format(input_sam)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input sam', dest='input')
    parser.add_argument('--output', help='output sam', dest='output')

    args = parser.parse_args()

    MarkDup(args.input, args.output)
