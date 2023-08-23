# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 11:45:25 2021

@author: FISXV1

This script takes an input path to a folder containing multiple synapse data session folders; all sessions within the folder must use the same
experimental conditions

THERE IS NO DOWNSAMPLING WITHIN THIS SCRIPT; CONSIDER ADDING WHEN NECESSARY
"""


#%% In [1] Import stuff

import os
import tdt
import numpy as np
import matplotlib.pyplot as plt

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
def load_sessions(path = r"/Users/shivasenthilkumar/Desktop/Robinson_Lab/Raw_Data/HighLight:LowLight"): 

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
    
    for d1 in range(0, len(data_sessions)):
        for d2 in range(0, len(data_sessions[d1].streams[sensor].filtered)):
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
epoch_ISOS_data_reformat = epoch_ISOS_data.reshape(np.shape(epoch_ISOS_data)[0],-1)
epoch_matched_data_reformat = e ed_data.reshape(np.shape(epoch_matched_data)[0], -1)


epoch_avg_streams, epoch_avg_ISOS = avg_across_sessions(epoch_ISOS_data_reformat, epoch_matched_data_reformat, dLight)
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
        
final_sig = np.empty((len(epoch_avg_stream)))
final_ISOS = np.empty((len(epoch_avg_streams)))

# subtract dc from both signals to normalize to zero, this places both signals ontop of one another
for n in range(0, len(epoch_avg_streams)):
    final_sig[n] = epoch_avg_streams[n]-dc_sig[n]
    final_ISOS[n] = epoch_avg_ISOS[n]-dc_ISOS[n]


#%% In [6] Set up some user defined variables for plot labeling - eventually, all labeling should be completely automatic

exp_con = input('What is your experimental condition?\n\n')
sessions = str(len(data_sessions))

    
#%% In [6] Plot the epoch averaged streams

t = int(input('which trial would you like to average across?\n\n')) - 1
fs = data_sessions[0].streams[dLight].fs

# Create the time vector for each stream store
ts1 = TRANGE[0] + np.linspace(1, len(final_sig[t]), len(final_sig[t]))/fs#*binsize - multiply fs by binsize if downsampling
ts2 = TRANGE[0] + np.linspace(1, len(final_ISOS[t]), len(final_ISOS[t]))/fs#*binsize - ^^^ same here

# Start making a figure with 4 subplots
# First plot is the 405 and 465 averaged signals
fig = plt.figure(figsize=(12, 22))
ax0 = fig.add_subplot(311) # work with axes and not current plot (plt.)

# Plotting the traces
p1, = ax0.plot(ts1, final_sig[t], linewidth=2, color='green', label="dLight")
p2, = ax0.plot(ts2, final_ISOS[t], linewidth=2, color='blueviolet', label="ISOS")

# Plotting standard error bands
p3 = ax0.fill_between(ts1, final_sig[t]+std_sig[t], final_sig[t]-std_sig[t],
                      facecolor='green', alpha=0.2)
p4 = ax0.fill_between(ts2, final_ISOS[t]+std_ISOS[t], final_ISOS[t]-std_ISOS[t],
                      facecolor='blueviolet', alpha=0.2)

# Plotting a line at t = 0
p5 = ax0.axvline(x=0, linewidth=3, color='slategray', label='Light Onset')

# Finish up the plot
ax0.set_xlabel('Seconds')
ax0.set_ylabel('mV')

# Use this title if adding artifact removal back in
# =============================================================================
# ax0.set_title(', %i Trials (%i Artifacts Removed)') #set this title to something that makes sense
#               % (len(data.streams[sensor].filtered), total_artifacts))
# =============================================================================

ax0.set_title(exp_con+ ' Stimulus; ' + sessions + ' Sessions in Data Set')

ax0.legend(handles=[p1, p2, p5], loc='upper right')
ax0.set_ylim(min(np.min(final_sig[t]-std_sig[t]), np.min(final_ISOS[t]-std_ISOS[t])),
             max(np.max(final_sig[t]+std_sig[t]), np.max(final_ISOS[t]+std_ISOS[t])))
ax0.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0]);

 # Jupyter cells will output any figure calls made, so if you don't want to see it just yet, close existing axis
            # https://stackoverflow.com/questions/18717877/prevent-plot-from-showing-in-jupyter-notebook
            # Note that this is not good code practice - Jupyter lends it self to these types of bad workarounds 
            
# =============================================================================
# #%% In [7]
# 
# def find_local_max(signal, time):
# 
#     #time = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]) #time vector array
#     #signal =  np.random.rand(15, 20000) #signal vector array
#     num_rows = np.shape(signal)[0] #number of rows for signal vector
#     Index_of_TTL = np.where(abs(time) == np.min(abs(time))) #for the 0, enter in the position of the x-coordinate of TTL signal
#     Index_of_TTL = int(Index_of_TTL[0]) #converts the 'tuple' of 1 element that is returned to an integer
#     Index_Final = np.where(abs(time-3)==np.min(abs(time-3))) #for the 5, enter in the position of the final time stamp
#     Index_Final = int(Index_Final[0]) #converts the 'tuple' of 1 element that is returned to an integer
#     rngArray = signal[0:num_rows, Index_of_TTL:Index_Final] #creates a new array (subarray) from the original array
#     localMax_percolumn = np.max(rngArray,axis=1) #find the maximum of each row in the subarray
#     trials = np.arange(1, (num_rows +1), 1) #creates the x-axis evenly spaced by 1 --> for the # of trials
# 
#     ax1 = fig.add_subplot(312)
#     ax1.set_ylim(0, 4)
#     ax1.set_xlim(0, 15)
#     scatterplot = ax1.scatter(trials, localMax_percolumn) #creates a scatter plot between trials and 
#   
#     return(scatterplot)
# #%%
# 
# 
# final_sig = np.array(final_sig)
# scatterplot = find_local_max(final_sig, ts1)
# =============================================================================

#%% Fit 405 channel onto 465 channel to detrend signal bleaching

Y_fit_all = []
Y_dF_all = []
for x, y in zip(epoch_ISOS_data_reformat, epoch_matched_data_reformat): #contain all 15 traces which are snips around the epoch data
    x = np.array(x)
    y = np.array(y)
    bls = np.polyfit(x, y, 1)
    fit_line = np.multiply(bls[0], x) + bls[1]
    Y_fit_all.append(fit_line)
    Y_dF_all.append(y-fit_line)

# Getting the z-score and standard error
zall = []
for dF in Y_dF_all: 
   ind = np.where((np.array(ts2)<BASELINE_PER[1]) & (np.array(ts2)>BASELINE_PER[0]))
   zb = np.mean(dF[ind])
   zsd = np.std(dF[ind])
   zall.append((dF - zb)/zsd)

zerror = np.std(zall, axis=0)/np.sqrt(np.size(zall, axis=0))

#%% Plot the z-score for 465 with std error bands

ax2 = fig.add_subplot(313)
x_axis = ts2
y_axis = np.mean(zall, axis=0)
max_zscore = y_axis.max(axis = 0)
min_zscore = y_axis.min(axis = 0)
p6 = ax2.plot(x_axis, y_axis, linewidth=2, color='green', label='dLight') #line graph
p7 = ax2.fill_between(x_axis, np.mean(zall, axis=0)+zerror
                       ,np.mean(zall, axis=0)-zerror, facecolor='green', alpha=0.2)
p8 = ax2.axvline(x=0, linewidth=3, color='slategray', label='Light Onset')
ax2.set_ylabel('z-Score')
ax2.set_xlabel('Seconds')
ax2.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0]) #user must adjust the min and max values to determine the range they want to see the bottom graph
#original is TRANGE[0], TRANGE[1]+TRANGE[0]
ax2.set_title('Light Stim Response')

#%% Heat map based on z-score of 405 fit subtracted 465

ax1 = fig.add_subplot(312)
cs = ax1.imshow(zall, cmap=plt.cm.plasma, interpolation='none', aspect="auto", vmin = (min_zscore - 0.75), vmax = (max_zscore + 0.75),
                extent=[TRANGE[0], TRANGE[1]+TRANGE[0], 0, len(data.streams[dLight].filtered)])
#user must adjust the min and max values to determine the range they want to see for the heat map graph
#original is TRANGE[0], TRANGE[1]+TRANGE[0]
cbar = fig.colorbar(cs, pad=0.01, fraction=0.02)

ax1.set_title('Individual z-Score Traces')
ax1.set_ylabel('Trials')
ax1.set_xlabel('Seconds from Shock Onset')

# Suppress figure output again
