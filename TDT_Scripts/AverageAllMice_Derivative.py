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
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
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
def load_sessions(path = r"/Users/shivasenthilkumar/Desktop/Robinson_Lab/Scripts/TDT_Scripts/DATA/High Light Intensity TTLS"): 

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


#%% In [6] Set up some user defined variables for plot labeling - eventually, all labeling should be completely automatic

exp_con = input('What is your experimental condition?\n\n')
sessions = str(len(data_sessions))

    
#%% In [7] Plot the epoch averaged streams


fs = data_sessions[0].streams[dLight].fs

# Create the time vector for each stream store
ts1 = TRANGE[0] + np.linspace(1, len(final_sig), len(final_sig))/fs#*binsize - multiply fs by binsize if downsampling
ts2 = TRANGE[0] + np.linspace(1, len(final_ISOS), len(final_ISOS))/fs#*binsize - ^^^ same here

# Start making a figure with 4 subplots
# First plot is the 405 and 465 averaged signals
fig = plt.figure(figsize=(20, 40))
# ax0 = fig.add_subplot(311) #work with axes and not current plot (plt.)

# # Plotting the traces
# p1, = ax0.plot(ts1, final_sig, linewidth=2, color='green', label="dLight")
# p2, = ax0.plot(ts2, final_ISOS, linewidth=2, color='blueviolet', label="ISOS")

# # Plotting standard error bands
# p3 = ax0.fill_between(ts1, final_sig+std_sig, final_sig-std_sig,
#                       facecolor='green', alpha=0.2)
# p4 = ax0.fill_between(ts2, final_ISOS+std_ISOS, final_ISOS-std_ISOS,
#                       facecolor='blueviolet', alpha=0.2)

# # Plotting a line at t = 0
# p5 = ax0.axvline(x=0, linewidth=3, color='slategray', label='Light Onset')

# # Finish up the plot
# ax0.set_xlabel('Seconds')
# ax0.set_ylabel('mV')

# # Use this title if adding artifact removal back in
# # =============================================================================
# # ax0.set_title(', %i Trials (%i Artifacts Removed)') #set this title to something that makes sense
# #               % (len(data.streams[sensor].filtered), total_artifacts))
# # =============================================================================

# ax0.set_title(exp_con+ ' Stimulus; ' + sessions + ' Sessions in Data Set')

# ax0.legend(handles=[p1, p2, p5], loc='upper right')
# ax0.set_ylim(min(np.min(final_sig-std_sig), np.min(final_ISOS-std_ISOS)),
#               max(np.max(final_sig+std_sig), np.max(final_ISOS+std_ISOS)))
# ax0.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0]);

 # Jupyter cells will output any figure calls made, so if you don't want to see it just yet, close existing axis
            # https://stackoverflow.com/questions/18717877/prevent-plot-from-showing-in-jupyter-notebook
            # Note that this is not good code practice - Jupyter lends it self to these types of bad workarounds 

#%% In [8] Z-Score + Find Standard Error Bars (Averaged Trace across All Mice)

matched_data = np.transpose (epoch_matched_data, (2, 1, 0))
matchedlist = matched_data.tolist()

zscore_all = []
zscoreavg_stream = []
zscorelocalmax = []
std_error = []

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
    
zscoreavg_stream = np.array(zscoreavg_stream)
zscoreavg_stream_new = zscoreavg_stream.transpose(1,0)

for x in zscoreavg_stream_new:
    stddev = np.std (x)
    ind_stderror = (stddev/(math.sqrt(np.shape(epoch_matched_data)[1])))
    std_error.append (ind_stderror)
    
avg_zscore_stream = np.mean(zscoreavg_stream_new, axis = 1)

#%% In [9] Find Time TTL

TTL_to_peak_total_time = []

for sessionrow in zscoreavg_stream:
    sessionrowmax = sessionrow.max(axis = 0)
    x = np.where(sessionrow == sessionrowmax)
    x = int(x[0])
    TTL_to_peak_ind = ts2[x]
    TTL_to_peak_total_time.append(TTL_to_peak_ind)
    
#%% In [10] Calculations needed for the heat map 

epoch_matched_data_rearrange = epoch_matched_data.transpose(1, 2, 0)
epoch_data_average = np.mean(epoch_matched_data_rearrange, axis = 2) #calculates the average of the row (per trials)
epoch_data_average = epoch_data_average.flatten()
epoch_data_std = np.std(epoch_matched_data_rearrange, axis = 2) #calculates the standard deviation of the row (per trials)
epoch_data_std = epoch_data_std.flatten()

zscore_list = []
modified_zscore_list = []

y = 0
for each_face in epoch_matched_data_rearrange:
    for each_row in each_face:
        zscore = ((each_row - epoch_data_average[y])/(epoch_data_std[y]))
        zscore_list.append(zscore)
        y += 1
        
for i in range (0, np.shape(epoch_matched_data)[2]):
    individuallist = zscore_list[i::np.shape(epoch_matched_data)[2]]
    modified_zscore_list.append(individuallist)
    
modified_zscore_list_array = np.array(modified_zscore_list)
new_array_reshape_rearrange_row_means = np.mean(modified_zscore_list_array, axis = 1)

#%% In [11] Plot the z-score for 465 with std error bands

ax2 = fig.add_subplot(511)
x_axis = ts2
y_axis = avg_zscore_stream
max_zscore = y_axis.max(axis = 0)
min_zscore = y_axis.min(axis = 0)
p6 = ax2.plot(x_axis, y_axis, linewidth=2, color='green', label='dLight') #line graph
p7 = ax2.fill_between(x_axis, avg_zscore_stream+std_error
                         ,avg_zscore_stream-std_error, facecolor='green', alpha=0.2)
p8 = ax2.axvline(x=0, linewidth=3, color='slategray', label='Light Onset')
ax2.set_ylabel('z-Score')
ax2.set_xlabel('Seconds')
ax2.set_xlim(TRANGE[0], TRANGE[1]+TRANGE[0]) #user must adjust the min and max values to determine the range they want to see the bottom graph
#original is TRANGE[0], TRANGE[1]+TRANGE[0]
ax2.set_title(exp_con)


#%% In [12] Heat map based on z-score of 405 fit subtracted 465

ax1 = fig.add_subplot(512)
cs = ax1.imshow(new_array_reshape_rearrange_row_means, cmap=plt.cm.plasma, interpolation='none', aspect="auto", vmin = min_zscore-0.5, vmax = max_zscore+0.5,
                extent=[TRANGE[0], TRANGE[1]+TRANGE[0], 0,  len(new_array_reshape_rearrange_row_means)])
#user must adjust the min and max values to determine the range they want to see for the heat map graph
#original is TRANGE[0], TRANGE[1]+TRANGE[0]
cbar = fig.colorbar(cs, pad=0.01, fraction=0.02)

ax1.set_title('Individual z-Score Traces')
ax1.set_ylabel('Trials')
ax1.set_xlabel('Seconds from Shock Onset')

#%% In [13] Derivative Calculation of Z-Score Trace

#this for loop goes through the average z scored stream and find the time Index of the TTL (i.e what index is the TTL located at)
#and finds the Final Index (i.e. what index is the desired final x (time) value). After the for loop, a new subarray of both the x 
#and y values from the TTL to a certain number of seconds is created.
for z in avg_zscore_stream:
    Index_of_TTL = np.where(abs(ts2) == np.min(abs(ts2))) 
    Index_of_TTL = int(Index_of_TTL[0]) 
    Index_Final = np.where(abs(ts2-1)==np.min(abs(ts2-1))) 
    Index_Final = int(Index_Final[0]) 
newtmArray = ts2[Index_of_TTL:(Index_Final+1)]
rngArray = avg_zscore_stream[Index_of_TTL:(Index_Final+1)] 

#finds the index of the maximum value of the function
Index_of_Max = np.where(rngArray == np.max(rngArray))
Index_of_Max = int(Index_of_Max[0])

#creates two new subarray from 0 (which is the TTL) to the time that corresponds to the Max Value
#this is true for both the x values (ttime) and y values (begtopeak)
begtopeak = rngArray[0:(Index_of_Max+1)]
ttime = newtmArray[0:(Index_of_Max+1)]

#finds all the indexes of the local minima, and then picks the maximumindex (because the maximum index of all the 
#minimum values is the beginning of the peak starting to spike)
minvalues_indexes = argrelextrema(begtopeak, np.less)
maxElement_index = np.amax(minvalues_indexes)
minvalue = begtopeak[maxElement_index]

#creates a horizontal line with the y value equalling the minvalue
horizonalline = minvalue*np.ones(len(newtmArray))
firstminindex = maxElement_index

#finds all the intersection points that the horizontal line passed through past the peak value 
#and gives back the second minimum value with lowest index value
outside_index_of_max = rngArray[Index_of_Max:-1]
subhorizontalline = horizonalline[Index_of_Max:-1]
idx = np.argwhere(np.diff(np.sign(subhorizontalline - outside_index_of_max))).flatten()

#calculates the secondminindex and the secondminvalue
secondminindex = Index_of_Max+idx[0] #based off TTL to Index Final (newtm and rngArray)
secondmin = rngArray[secondminindex] #returns the second minimum value  

#creates the final time and zscore functions to be plotted
finaltm = newtmArray[0:secondminindex]
finaltm = 1000*finaltm
finalyvalue = rngArray[0:secondminindex]

#calculated the 1st and 2nd Derivative of the original function
derivative = np.gradient(finalyvalue)
secondderivative = np.gradient(derivative)

#creating the new figure for the graphs to be plotted on
fig1, host = plt.subplots(figsize=(12,15))
par1 = host.twinx()
par2 = host.twinx()
    
#plotting on same graph
host.set_xlabel("Time (ms)") 
host.set_xlim(0, finaltm[-1])
# host.set_ylabel("Average Z-Score Trace")
host.set_ylim(-(np.max(finalyvalue)+0.2), np.max(finalyvalue)+0.2)
# par1.set_ylabel("1st Derivative of Avg Z-Score Trace (zscore/s)")
par1.set_ylim(-(np.max(derivative)+0.02), np.max(derivative)+0.02)
# par2.set_ylabel("2nd Derivative of Avg Z-Score Trace (zscore/s^2)")
par2.set_ylim((np.min(secondderivative)-0.001), -(np.min(secondderivative)-0.001))

#the colors for each individual lien
color1 = plt.cm.viridis(0.2670)
color2 = plt.cm.viridis(0.4264)
color3 = plt.cm.viridis(0.7812)

#plotting each individual function
p1, = host.plot(finaltm, finalyvalue, linewidth=2, color=color1, label="Average Z-Score Trace")
p2, = par1.plot(finaltm, derivative, linewidth=2, color=color2, label="1st Derivative (zscore/sec)")
p3, = par2.plot(finaltm, secondderivative, linewidth=2, color=color3, label="2nd Derivative (zscore/sec^2)")
p4 = par1.axhline(y=0, color='r', linestyle='--')

lns = [p1, p2, p3]
host.legend(handles=lns, loc='best')

# no y-ticks             
host.yaxis.set_ticks([])    
par1.yaxis.set_ticks([])    
par2.yaxis.set_ticks([])

#the spacing of the ticks on the x-axis
tick_spacing = 50
host.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())

fig.tight_layout()

#calculates the max-peak values and their corresponding time_values
peak1 = np.max(finalyvalue)
peak2 = np.max(derivative)
peak3 = np.max(secondderivative)

peak1_index = int (np.where(finalyvalue == peak1)[0])
peak2_index = int (np.where(derivative == peak2)[0])
peak3_index = int (np.where(secondderivative == peak3)[0])

peak1_time = finaltm[peak1_index]
peak2_time = finaltm[peak2_index]
peak3_time = finaltm[peak3_index]



