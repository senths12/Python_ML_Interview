# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 11:45:25 2021

@author: Shiva

This script takes an input path to a folder containing multiple synapse data session folders; all sessions within the folder must use the same
experimental conditions

THERE IS NO DOWNSAMPLING WITHIN THIS SCRIPT; CONSIDER ADDING WHEN NECESSARY
"""


#%% In [1] Import stuff

import os
import tdt
import numpy as np

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
        print (paths[x])
            
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
    
    for d1 in range(0, len(data_sessions)):
        for d2 in range(0, (len(data_sessions[d1].streams[sensor].filtered))):
           epoch_matched_data[:, d1, d2] = data_sessions[d1].streams[sensor].filtered[d2]
           epoch_ISOS_data[:, d1, d2] = data_sessions[d1].streams[ISOS].filtered[d2]
    
            
    return((epoch_matched_data, epoch_ISOS_data))

# Create a list of epoch streams averaged across sessions by trial number; this takes the epoch_matched_data 3D array and processes down 
def avg_across_sessions(epoch_matched_data, epoch_ISOS_data, sensor, ISOS=ISOS):
    
    epoch_avg_streams = epoch_matched_data.mean(axis=1)
    epoch_avg_ISOS = epoch_ISOS_data.mean(axis=1)
    
    
    
    return(epoch_avg_streams, epoch_avg_ISOS)


#%% In [5] This calls the functions and ties them together

# Load all data sessions within the path and create the epoch streams within the sessions data structures
data_sessions = load_sessions()

for d in range(0, len(data_sessions)):
    data_sessions[d] = tdt.epoc_filter(data_sessions[d], REF_EPOC, t=TRANGE, values=TTL_CODE)

# match the length of all epoch streams - streams can vary in size by 1; 
# this finds the smaller size and clips off an extra value from the larger snips; returns the dLight/GCaMP snip size and the ISOS snip size
min1, min2 = epoch_size_match(data_sessions, dLight)


epoch_matched_data, epoch_ISOS_data = epoch_match(data_sessions, dLight)


#reshape the 3D collection of all data streams (epoch_matched_data) into a 2D array to average across all epoch filtered streams (trials and sessions/mice)
epoch_ISOS_data_reformat = epoch_ISOS_data.reshape(np.shape(epoch_ISOS_data)[0],-1)
epoch_matched_data_reformat = epoch_matched_data.reshape(np.shape(epoch_matched_data)[0], -1)



epoch_avg_streams, epoch_avg_ISOS = avg_across_sessions(epoch_matched_data_reformat, epoch_ISOS_data_reformat, dLight)
epoch_avg_streams = epoch_avg_streams.transpose()
epoch_avg_ISOS = epoch_avg_ISOS.transpose()


dc_sig = []
dc_ISOS = []
std_sig = []
std_ISOS = []


for n in range(0, len(epoch_avg_streams)):
    dc_sig.append(np.mean(epoch_avg_streams[n:]))
    dc_ISOS.append(np.mean(epoch_avg_ISOS[n:]))
    std_sig.append(np.std(epoch_avg_streams[n:])/np.sqrt(len(epoch_avg_streams[n:])))
    std_ISOS.append(np.std(epoch_avg_ISOS[n:])/np.sqrt(len(epoch_avg_ISOS[n:])))
        
final_sig = np.empty((len(epoch_avg_streams)))
final_ISOS = np.empty((len(epoch_avg_streams)))

# subtract dc from both signals to normalize to zero, this places both signals ontop of one another
for n in range(0, len(epoch_avg_streams)):
    final_sig[n] = epoch_avg_streams[n]-dc_sig[n]
    final_ISOS[n] = epoch_avg_ISOS[n]-dc_ISOS[n]

#%% In [6] Plot the epoch averaged streams


fs = data_sessions[0].streams[dLight].fs

# Create the time vector for each stream store
ts1 = TRANGE[0] + np.linspace(1, len(final_sig), len(final_sig))/fs#*binsize - multiply fs by binsize if downsampling
ts2 = TRANGE[0] + np.linspace(1, len(final_ISOS), len(final_ISOS))/fs#*binsize - ^^^ same here

matched_data = np.transpose (epoch_matched_data, (2, 1, 0))
matchedlist = matched_data.tolist()

zscore_all = []
zscoreavg_stream = []
zscorelocalmax = []


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
    Index_of_TTL = np.where(abs(ts2) == np.min(abs(ts2))) 
    Index_of_TTL = int(Index_of_TTL[0]) 
    Index_Final = np.where(abs(ts2-2)==np.min(abs(ts2-2))) 
    Index_Final = int(Index_Final[0]) 
    rngArray = z[Index_of_TTL:Index_Final] 
    zscorelocalmax.append(np.max(rngArray))


