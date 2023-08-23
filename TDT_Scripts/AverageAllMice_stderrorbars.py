# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:32:49 2021

@author: auste
"""

#%% In [1] Import necessary packages

import os
import tdt
import numpy as np
import matplotlib.pyplot as plt
import math

#%% In [2] plt imported separately due to Jupyter bug

import matplotlib
matplotlib.rcParams['font.size'] = 16 #set font size for all plots

#%% In [3] Set up variables

REF_EPOC = 'PtC0' #where the ttl data is stored; holds timestamps for ttl signals

TTL_CODE = [1] #shock onset event code we are interested in; this is the onset time of the ttl

# make some variables up here to so if they change in new recordings you won't
# have to change everything downstream
ISOS = '_405A' # 405nm channel. Formally STREAM_STORE1 in maltab example; change to _405A for our system
dLight = '_465A' # 465nm channel. Formally STREAM_STORE2 in maltab example; change to _465A for our system
TRANGE = [-2, 16] # window size [start time relative to epoc onset, window duration] # works when set to [-10, 20], small adjustments work but breaks with larger changes (i.e. [-10, 21 through 24] works, [-10, 25] suddenly doesn't)
# this may be an issue where the data stream cuts off shortly after the final light delivery - probably more issues than this
BASELINE_PER = [-10, -1] # baseline period within our window; reference period for calculating the z-score
ARTIFACT = float("inf") # optionally set an artifact rejection level

#call read block - new variable 'data' is the full data structure

#%% In [4] Define functions for data manipulation

# Load all data folders within a parent folder into list called data_sessions
def load_sessions(path = r"/Users/shivasenthilkumar/Desktop/Robinson_Lab/Scripts/TDT_Scripts/DATA/Darkness TTLs"): 

    folders = os.listdir(path)
    folders.remove('.DS_Store')
    
    paths = []
    data_sessions = []
    
    
    #create a list of strings of the session folder paths
    for x in range(0, len(folders)):
        p = (path + "//" + folders[x]) #This plus the following line can be easily condensed into one
       # p = p.replace("\\", "/")
        paths.append(p)
            
    #create a list of all session data structures for analysis
    for x in range(0, len(paths)):
        d = tdt.read_block(paths[x])
        data_sessions.append(d)
        
    return(data_sessions)


# Match the size of all epoch snips from all sessions   
def epoch_size_match(data_sessions, sensor, ISOS=ISOS): #sensor is a string pointing to the desired data stream (usually GCaMP/dLight = _465A)

    epoch_mins = []
    isos_mins = []

    for data in data_sessions:
        min1 = np.min([np.size(x) for x in data.streams[sensor].filtered])
        min2 = np.min([np.size(x) for x in data.streams[ISOS].filtered]) 
        epoch_mins.append(min1)
        isos_mins.append(min2)
        
    for data in data_sessions:
        data.streams[sensor].filtered = [x[1:min1] for x in data.streams[sensor].filtered]
        data.streams[ISOS].filtered = [x[1:min2] for x in data.streams[ISOS].filtered]
        
    return([min1, min2])

# Downsample to desired bin size (default setting is 10); this needs to be modified for averaging across sessions; this should be done within
# the data structure themselves, rather than compiled into python lists as seen here
# =============================================================================
# def downsample(data_sessions, sensor, min1, min2, ISOS = ISOS, binsize = 10):
#     
#     for data in data_sessions:
#         
#         N = binsize # Average every 10 samples into 1 value
#         F405 = []
#         F465 = []
#         for lst in data.streams[ISOS].filtered: 
#             small_lst = []
#             for i in range(0, min2, N):
#                 small_lst.append(np.mean(lst[i:i+N-1])) # This is the moving window mean
#             F405.append(small_lst)
#         
#         for lst in data.streams[sensor].filtered: 
#             small_lst = []
#             for i in range(0, min1, N):
#                 small_lst.append(np.mean(lst[i:i+N-1]))
#             F465.append(small_lst)
#         
#     return([F405, F465, binsize])
# =============================================================================

# Collect data in format for averaging across sessions by trial number; creates a 3D array called epoch_matched_data; rows within each 2D slice
# of the array contain individual snips within each row; all snips within the 2D subarray are from the same trial #, but different sessions
def epoch_match(data_sessions, sensor, ISOS = ISOS):
    
    #Initiate some np.arrays to store data in
    
    epoch_matched_data = np.empty((len(data_sessions[0].streams[sensor].filtered[0]), len(data_sessions), len(data_sessions[0].streams[sensor].filtered)))
    epoch_ISOS_data = np.empty((len(data_sessions[0].streams[sensor].filtered[0]), len(data_sessions), len(data_sessions[0].streams[sensor].filtered)))
#
#
# GET RID OF THE -1 IN THE FOLLOWING LOOP, second line -1 added to account for missing last epoch
#
#
    
    for d1 in range(0, len(data_sessions)):
        for d2 in range(0, len(data_sessions[d1].streams[sensor].filtered)): #subtract 1 from the length for the tone data (at least one of the streams must be missing the list epoch)
            epoch_matched_data[:, d1, d2] = data_sessions[d1].streams[sensor].filtered[d2]
            epoch_ISOS_data[:, d1, d2] = data_sessions[d1].streams[ISOS].filtered[d2]
    
#
#
#
#
#            
    return((epoch_matched_data, epoch_ISOS_data))

# Create an array of epoch streams averaged across trials within each session; this takes the epoch_matched_data 3D array and processes down 
def avg_across_sessions(epoch_matched_data, epoch_ISOS_data): #sensor, and ISOS=ISOS removed from inputs, these should not be needed after the epoch_size_match function
    
    epoch_avg_streams = np.mean(epoch_matched_data, axis=2) 
    epoch_avg_ISOS = np.mean(epoch_matched_data, axis=2) 
    
    
    
    return(epoch_avg_streams, epoch_avg_ISOS)

# z-score all streams within a 2D numpy array in form (len stream x # of sessions)
#def z-streams(streams_array):    #put the z-scoring within sessions function here when complete
    
    



#%% In [5] This calls the functions and ties them together

# Load all data sessions within the path and create the epoch streams within the sessions data structures
data_sessions = load_sessions()

for d in range(0, len(data_sessions)):
    data_sessions[d] = tdt.epoc_filter(data_sessions[d], REF_EPOC, t=TRANGE, values=TTL_CODE)

# match the length of all epoch streams - streams can vary in size by 1; 
# this finds the smaller size and clips off an extra value from the larger snips; returns the dLight/GCaMP snip size and the ISOS snip size
min1, min2 = epoch_size_match(data_sessions, dLight)

#Get all raw data streams and ISOS streams into a 3D numpy array (organized by trial and session)
epoch_matched_data, epoch_ISOS_data = epoch_match(data_sessions, dLight)

#Average the stream and ISOS 3D arrays across trials - returns a 2D array of the epoch filtered and session averaged data
epoch_avg_sessions, epoch_avg_ISOS = avg_across_sessions(epoch_matched_data, epoch_ISOS_data)

#%% In [6] Z-Score + Find Standard Error Bars

epoch_matched_data_rearrange = epoch_matched_data.transpose(1, 2, 0)
avgstream_session = np.mean(epoch_matched_data_rearrange, axis = 1)

zscore_zall = []
zsd_list = []
std_error = []
for eachsession in avgstream_session: 
    zb = np.mean (eachsession)
    zsd = np.std (eachsession)
    zscore_zall.append((eachsession - zb)/zsd)
    zsd_list.append(zsd)
for x in zsd_list:
    ind_stderror = (x/(math.sqrt(np.shape(epoch_matched_data)[2])))
    std_error.append(ind_stderror)
zscore_zall = np.array(zscore_zall)
avg_zscore_stream = np.mean(zscore_zall, axis = 0)
zerror = np.mean(std_error, axis = 0)


    











