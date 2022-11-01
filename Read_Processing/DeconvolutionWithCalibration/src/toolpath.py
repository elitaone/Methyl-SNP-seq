# -*- coding: utf-8 -*-
# python2.7
'''
Test samtools and bedtools
'''

import os
from subprocess import check_output

try:
    import pandas
    import pathos
except:
    print 'Need pandas and pathos module installed.'
    quit()

def init(samtools, bedtools, bowtie2):
    '''set tool path as global variables'''
    global samtoolspath
    samtoolspath = samtools
    global bedtoolspath
    bedtoolspath = bedtools
    global bowtie2path
    bowtie2path = bowtie2

def Tools():
    '''
    Test whether bedtools or samtools works
    '''
    with open('mytools.version', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{} --version'.format(bedtoolspath)
    try:
        output1 = check_output(['sh', 'mytools.version'], shell=False)
        os.remove('mytools.version')
        print "Find bedtools"
    except:
        print 'bedtools error.'
        os.remove('mytools.version')
        quit()
    
    with open('mytools.version', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{} --version'.format(samtoolspath)
    try:
        output2 = check_output(['sh', 'mytools.version'], shell=False)
        os.remove('mytools.version')
        print "Find samtools"
    except:
        print 'samtools error.'
        os.remove('mytools.version')
        quit()

    with open('mytools.version', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{} --version'.format(bowtie2path)
    try:
        output3 = check_output(['sh', 'mytools.version'], shell=False)
        os.remove('mytools.version')
        print "Find bowtie2"
    except:
        print 'bowtie2 error.'
        os.remove('mytools.version')
        quit()
    
    return output1, output2, output3
