# -*- coding: utf-8 -*-
"""
Created on Web Jun 30 

@author: shiva
"""

#%% In [1] Import necessary packages

import os
import tdt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import xlwt
from xlwt import Workbook


mmatched_data = np.transpose (epoch_matched_data, (2, 1, 0))
matchedlist = matched_data.tolist()

zscore_all = []
zscoreavg_stream = []
zscorelocalmax = []
zscore_stderror = []

#finds all the zscores
for triallist in matchedlist: #for mouse n in trial n
    for sublist in triallist: #for each data stream in mouse n
         zb = np.mean (sublist)
         zsd = np.std (sublist)
         zscore_all.append((sublist - zb)/zsd)
zscore_all = np.array(zscore_all)
total = np.shape(epoch_matched_data)[1]*np.shape(epoch_matched_data)[2]
length = np.shape(epoch_matched_data)[1]
width = np.shape(epoch_matched_data)[2]
original_array = zscore_all[0:total:length]
i = 0
for i in range(1, length):
    new_array = zscore_all[i:total:length]
    total_array = np.concatenate((original_array, new_array))
    original_array = total_array
    i += 1
for y in range(0, total, width):
    zscoreavg_stream.append(np.mean((total_array[y:(y+width)-1]), axis = 0))
for z in zscoreavg_stream:
    Index_of_TTL = np.where(abs(ts1) == np.min(abs(ts1))) 
    Index_of_TTL = int(Index_of_TTL[0]) 
    Index_Final = np.where(abs(ts1-2)==np.min(abs(ts1-2))) 
    Index_Final = int(Index_Final[0]) 
    rngArray = z[Index_of_TTL:Index_Final] 
    zscorelocalmax.append(np.max(rngArray))
         

