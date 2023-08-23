# -*- coding: utf-8 -*-
"""
Created on Fri May 14 14:20:30 2021

@author: FISXV1
"""

#%%
import numpy as np
#from sklearn.metrics import auc #used in the commented cells (AUC and barplot); uneeded for now
import matplotlib.pyplot as plt  # standard Python plotting library
#import scipy.stats as stats #used in the commented cells; uneeded for now

import tdt

folderpath = ""
def readData (folderpath):
    fullstr = folderpath 
    substr = "Desktop"
    origsubstrlgnth = len (substr)
    partialstr = ""
    i = 0
    d_value = fullstr.find(substr)
    if substr in fullstr:
        index_of_data = origsubstrlgnth + d_value
        while (i + origsubstrlgnth < (len(fullstr)) - d_value):
            partialstr = substr + fullstr[index_of_data + i]
            substr = partialstr
            i += 1
    return (substr)

#%%
#BLOCKPATH = 'data\\24L-210512-105617'

BLOCKPATH = readData(input("Please enter the folder path name"))
data = tdt.read_block(BLOCKPATH)
title = data.info.blockname
