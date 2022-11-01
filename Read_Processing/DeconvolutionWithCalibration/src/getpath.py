# -*- coding: utf-8 -*-
# python2.7
'''
Contain useful functions
'''
import os
import time

def GetFilePath(input_file):
    '''
    return the abspath of the input_file
    if use relative path and input_file is not in cwd, the abspath will be wrong, return False
    '''
    filepath = os.path.abspath(input_file)
    
    if os.path.exists(filepath):
        return filepath
    else:
        print 'Can not find file: {}'.format(filepath)
        return None
        
def CreatePrefix():
    '''
    create prefix based on current time
    '''
    localtime = time.asctime(time.localtime())
    prefix = ''.join(localtime.split()[-2].split(':')) # '151542'
    return prefix

def CheckFormat(fileformat):
    '''
    check the format of fileformat in python2, compressed or unconmpressed
    '''
    GZIP_MAGIC_NUMBER = "1f8b"
    with open(fileformat) as f:
        if f.read(2).encode("hex") == GZIP_MAGIC_NUMBER:
            return True # gz
        else:
            return False # uncompressed