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
#%%
#BLOCKPATH = 'data\\24L-210512-105617'

BLOCKPATH = '/Users/shivasenthilkumar/Desktop/Robinson_Lab/Scripts/TDT_Scripts/DATA/Looming/148N-210810-NewContrastLoom.5Pulses.5Trials'
data = tdt.read_block(BLOCKPATH)
title = data.info.blockname

#%%
import matplotlib 
matplotlib.rcParams['font.size'] = 16 #set font size for all plots

#%% This is for assigning a unique code to every epoch event 

# for i in range(0, len(data.epocs.PtC0.data)):
   # data.epocs.PtC0.data[i] = i+1


#%% Set up variables

REF_EPOC = 'PC0/' #where the ttl data is stored; holds timestamps for ttl signals
TTL_CODE = [1] #shock onset event code we are interested in; this is the onset time of the ttl

# make some variables up here to so if they change in new recordings you won't
# have to change everything downstream
ISOS = '_405A' # 405nm channel. Formally STREAM_STORE1 in maltab example; change to _405A for our system
dLight = '_465A' # 465nm channel. Formally STREAM_STORE2 in maltab example; change to _465A for our system
TRANGE = [-4, 20] # window size [start time relative to epoc onset, window duration] # works when set to [-10, 20], small adjustments work but breaks with larger changes (i.e. [-10, 21 through 24] works, [-10, 25] suddenly doesn't)
# this may be an issue where the data stream cuts off shortly after the final light delivery - probably more issues than this
BASELINE_PER = [-10, 12] # baseline period within our window; reference period for calculating the z-score

#what is the baseline period window? standard deviation from where it is calculated from 

ARTIFACT = float("inf") # optionally set an artifact rejection level

#call read block - new variable 'data' is the full data structure


#%% TDT epoc_filter function to extract data around epoch event

data = tdt.epoc_filter(data, REF_EPOC, t=TRANGE, values=TTL_CODE) #this is taking in the different parameters

# Optionally remove artifacts. If any waveform is above ARTIFACT level, or
# below -ARTIFACT level, remove it from the data set.
total1 = np.size(data.streams[dLight].filtered)
total2 = np.size(data.streams[ISOS].filtered)

# List comprehension checking if any single array in 2D filtered array is > Artifact or < -Artifact
data.streams[dLight].filtered = [x for x in data.streams[dLight].filtered 
                                if not np.any(x > ARTIFACT) or np.any(x < -ARTIFACT)]

#whis is this saying?
data.streams[ISOS].filtered = [x for x in data.streams[ISOS].filtered 
                               if not np.any(x > ARTIFACT) or np.any(x < -ARTIFACT)]

#what is this saying?
# Get the total number of rejected arrays --> why do we need rejected arrays
bad1 = total1 - np.size(data.streams[dLight].filtered)
bad2 = total2 - np.size(data.streams[ISOS].filtered)
total_artifacts = bad1 + bad2

#%% Data segments may vary in length by 1 when applying a time filter; find the minimum and trim extra point in longer segments
#what is a minimum and trim extra point in longer segments?

# More examples of list comprehensions
min1 = np.min([np.size(x) for x in data.streams[dLight].filtered])
min2 = np.min([np.size(x) for x in data.streams[ISOS].filtered])
data.streams[dLight].filtered = [x[1:min1] for x in data.streams[dLight].filtered]

#what does the above code mean?

data.streams[ISOS].filtered = [x[1:min2] for x in data.streams[ISOS].filtered]

# Downsample and average 10x via a moving window mean
N = 10 # Average every 10 samples into 1 value
F405 = []
F465 = []
for lst in data.streams[ISOS].filtered: 
    small_lst = []
    for i in range(0, min2, N):
        small_lst.append(np.mean(lst[i:i+N-1])) # This is the moving window mean
    F405.append(small_lst)

for lst in data.streams[dLight].filtered: 
    small_lst = []
    for i in range(0, min1, N):
        small_lst.append(np.mean(lst[i:i+N-1]))
    F465.append(small_lst)

#Create a mean signal, standard error of signal, and DC offset (what is a DC offset)
meanF405 = np.mean(F405, axis=0)
stdF405 = np.std(F405, axis=0)/np.sqrt(len(data.streams[ISOS].filtered))
dcF405 = np.mean(meanF405)
meanF465 = np.mean(F465, axis=0)
stdF465 = np.std(F465, axis=0)/np.sqrt(len(data.streams[dLight].filtered))
dcF465 = np.mean(meanF465)

#%% Plot epoch averaged response --> gives my mV vs seconds response

# Create the time vector for each stream store
ts1 = TRANGE[0] + np.linspace(1, len(meanF465), len(meanF465))/data.streams[dLight].fs*N
#ts1 are all the seconds in the x-axis corresponding to the dLigth line

ts2 = TRANGE[0] + np.linspace(1, len(meanF405), len(meanF405))/data.streams[ISOS].fs*N
#ts2 are all the seconds in the x-axis corresponding to the ISOS line

# Subtract DC offset to get signals on top of one another
meanF405 = meanF405 - dcF405
meanF465 = meanF465 - dcF465

# Start making a figure with 4 subplots --> what is this subplot meaning?
# First plot is the 405 and 465 averaged signals


#normalize the 465 to the 405 --> what does that mean? for photobleaching
#405 is the control wavelength
#465 is the dLight signal
fig = plt.figure(figsize=(12, 22))
ax0 = fig.add_subplot(311) # work with axes and not current plot (plt.)

#where do you get the 311 from?

# Plotting the traces --> oh yes
p1, = ax0.plot(ts1, meanF465, linewidth=2, color='green', label='dlight')
p2, = ax0.plot(ts2, meanF405, linewidth=2, color='blueviolet', label='ISOS')

# Plotting standard error bands --> for the first graph
p3 = ax0.fill_between(ts1, meanF465+stdF465, meanF465-stdF465,
                      facecolor='green', alpha=0.2)
p4 = ax0.fill_between(ts2, meanF405+stdF405, meanF405-stdF405,
                      facecolor='blueviolet', alpha=0.2)

# Plotting a line at t = 0 --> plotting your ttl
p5 = ax0.axvline(x=0, linewidth=3, color='slategray', label='Light Onset')

# Finish up the plot --> what is total artifacts??
ax0.set_xlabel('Seconds')
ax0.set_ylabel('mV')
ax0.set_title(title + ', %i Trials (%i Artifacts Removed)'
              % (len(data.streams[dLight].filtered), total_artifacts))
ax0.legend(handles=[p1, p2, p5], loc='upper right')
ax0.set_ylim(min(np.min(meanF465-stdF465), np.min(meanF405-stdF405)),
             max(np.max(meanF465+stdF465), np.max(meanF405+stdF405)))

#what does this piece of code do? -->sets the zoom-in version for the first graph
ax0.set_xlim(TRANGE[0], TRANGE[0] + TRANGE[1]); #user must adjust the min and max values to determine the range they want to see in the graph
#original is TRANGE[0], TRANGE[1]

 # Jupyter cells will output any figure calls made, so if you don't want to see it just yet, close existing axis
            # https://stackoverflow.com/questions/18717877/prevent-plot-from-showing-in-jupyter-notebook
            # Note that this is not good code practice - Jupyter lends it self to these types of bad workarounds 
            
#%% Fit 405 channel onto 465 channel to detrend signal bleaching

Y_fit_all = []
Y_dF_all = []
for x, y in zip(F405, F465): #contain all 15 traces which are snips around the epoch data
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

#%% quantify changes as an area under the curve for cue (-5sec) vs shock (0sec)

# =============================================================================
# AUC = [] # cue, shock
# ind1 = np.where((np.array(ts2)<-3) & (np.array(ts2)>-5))
# AUC1= auc(ts2[ind1], np.mean(zall, axis=0)[ind1])
# ind2 = np.where((np.array(ts2)>0) & (np.array(ts2)<2))
# AUC2= auc(ts2[ind2], np.mean(zall, axis=0)[ind2])
# AUC.append(AUC1)
# AUC.append(AUC2)
# 
# # Run a two-sample T-test
# t_stat,p_val = stats.ttest_ind(np.mean(zall, axis=0)[ind1],
#                                np.mean(zall, axis=0)[ind2], equal_var=False)
# 
# =============================================================================

#%%

# =============================================================================
# ax3 = fig.add_subplot(414)
# p9 = ax3.bar(np.arange(len(AUC)), AUC, color=[.8, .8, .8], align='center', alpha=0.5)
# 
# # statistical annotation
# x1, x2 = 0, 1 # columns indices for labels
# y, h, col = max(AUC) + 2, 2, 'k'
# ax3.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# p10 = ax3.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)
# 
# # Finish up the plot
# ax3.set_ylim(0,y+2*h)
# ax3.set_ylabel('AUC')
# ax3.set_title('Cue vs Shock Response Changes')
# ax3.set_xticks(np.arange(-1, len(AUC)+1))
# ax3.set_xticklabels(['','Cue','Shock',''])
# 
# fig.tight_layout()
# 
# =============================================================================
fig


