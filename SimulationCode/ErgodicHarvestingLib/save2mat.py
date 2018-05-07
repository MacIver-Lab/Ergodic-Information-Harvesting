# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:12:32 2017

@author: Chen Chen
"""
# Path and file
from os import listdir
import fnmatch
# Import file IO modules
from scipy.io import savemat

def save2mat(fname, datDict, path='', appendIdx=True):
    if appendIdx:
        # Count files with the same name
        if len(path) > 0:
            fileCnt = len(fnmatch.filter(listdir(path), fname+'*'))
        else:
            fileCnt = len(fnmatch.filter(listdir('.'), fname+'*'))
        
        # Append file name with file count as an unique identifier
        savemat(path + fname + '-' + str(fileCnt), datDict)
    else:
        savemat(path + fname, datDict)
    
def checkMatFileExist(fname, path=''):
    if len(path) > 0:
        fileCnt = len(fnmatch.filter(listdir(path), fname+'*'))
    else:
        fileCnt = len(fnmatch.filter(listdir('.'), fname+'*'))
    
    return fileCnt