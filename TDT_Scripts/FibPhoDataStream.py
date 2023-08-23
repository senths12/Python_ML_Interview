# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:34:44 2021

@author: FISXV1
"""

#%% In [1]

import numpy as np
import matplotlib.pyplot as plt  # standard Python plotting library
#from io import StringIO

# import the tdt library
import tdt

#%% In [2]

BLOCKPATH =  "/Users/shivasenthilkumar/Desktop/Robinson_Lab/Scripts/TDT_Scripts/DATA/Other Loom/148N-210712-LargeLeftLoomingStim.3Loom.5Pulses"
data = tdt.read_block(BLOCKPATH)

#%% In [3]

import matplotlib
matplotlib.rcParams['font.size'] = 18 # set font size for all figures

# Make some variables up here to so if they change in new recordings you won't have to change everything downstream
GCAMP = '_465A' # GCaMP channel
ISOS = '_405A' # Isosbestic channel
LIGHT = 'PtC0'

#%% In [4] Plot raw demodulated response 

time = np.linspace(1,len(data.streams[GCAMP].data), len(data.streams[GCAMP].data))/data.streams[GCAMP].fs

# Plot both unprocessed demodulated stream            
fig1 = plt.figure(figsize=(10, 6))
ax0 = fig1.add_subplot(111)

# Plotting the traces
p1, = ax0.plot(time, data.streams[GCAMP].data, linewidth=2, color='green', label='dLight')
p2, = ax0.plot(time, data.streams[ISOS].data, linewidth=2, color='blueviolet', label='ISOS')

ax0.set_ylabel('mV')
ax0.set_xlabel('Seconds')
ax0.set_title('Raw Demodulated Response')
ax0.legend(handles=[p1,p2], loc='upper right')
fig1.tight_layout()

#%% In [5] Detrending the signal

# There is often a large artifact on the onset of LEDs turning on
# Remove data below a set time t
# t is the starting point of the figure
t = 30
inds = np.where(time>t)
ind = inds[0][0]

#set t_end to large number to capture all data at the end
t_end = 2500 #This allows for clipping from the end for when photometry was recorded beyond the intended end point
inds_end = np.where(time<t_end)
ind_end = inds_end[0][-1]


time = time[ind:ind_end] # go from ind to final index
data.streams[GCAMP].data = data.streams[GCAMP].data[ind:ind_end]
data.streams[ISOS].data = data.streams[ISOS].data[ind:ind_end]

# Plot again at new time range
fig2 = plt.figure(figsize=(10, 6))
ax1 = fig2.add_subplot(111)

# Plotting the traces
p1, = ax1.plot(time,data.streams[GCAMP].data, linewidth=2, color='green', label='dLight')
p2, = ax1.plot(time,data.streams[ISOS].data, linewidth=2, color='blueviolet', label='ISOS')

ax1.set_ylabel('mV')
ax1.set_xlabel('Seconds')
ax1.set_title('Raw Demodulated Response with Artifact Removed')
ax1.legend(handles=[p1,p2],loc='upper right')
fig2.tight_layout()

#%% In [6] Downsample, local averaging

N = 20 # Average every N samples into 1 value; poly fit poorly conditioned when this is low (e.g. 2-10 had issues)
F405 = []
F465 = []

for i in range(0, len(data.streams[GCAMP].data), N):
    F465.append(np.mean(data.streams[GCAMP].data[i:i+N-1])) # This is the moving window mean
data.streams[GCAMP].data = F465

for i in range(0, len(data.streams[ISOS].data), N):
    F405.append(np.mean(data.streams[ISOS].data[i:i+N-1]))
data.streams[ISOS].data = F405

#decimate time array to match length of demodulated stream
time = time[::N] # go from beginning to end of array in steps on N
time = time[:len(data.streams[GCAMP].data)]

# Detrending and dFF
# Full trace dFF according to Lerner et al. 2015
# http://dx.doi.org/10.1016/j.cell.2015.07.014
# dFF using 405 fit as baseline

x = np.array(data.streams[ISOS].data)
y = np.array(data.streams[GCAMP].data)
bls = np.polyfit(x, y, 1)
Y_fit_all = np.multiply(bls[0], x) + bls[1]
Y_dF_all = y - Y_fit_all


dFF = np.multiply(100, np.divide(Y_dF_all, Y_fit_all))
std_dFF = np.std(dFF)


#%%
# First make a continous time series of Licking TTL events (epocs) and plot
LIGHT_on = data.epocs[LIGHT].onset
LIGHT_off = []

for i in LIGHT_on:
    i = i + 10
    LIGHT_off.append(i)


# Add the first and last time stamps to make tails on the TTL stream
LIGHT_x = np.append(np.append(time[0], np.reshape(np.kron([LIGHT_on, LIGHT_off], 
                                                          np.array([[1], [1]])).T, [1,-1])[0]), time[-1])
sz = len(LIGHT_on)
d = data.epocs.PtC0.data
# Add zeros to beginning and end of 0,1 value array to match len of LICK_x
LIGHT_y = np.append(np.append(0,np.reshape(np.vstack([np.zeros(sz),
    d, d, np.zeros(sz)]).T, [1, -1])[0]),0)

y_scale = 10 #adjust according to data needs
y_shift = -20 #scale and shift are just for asthetics

#%% In [7] Plot dF/F

fig3 = plt.figure(figsize=(10,10))
ax2 = fig3.add_subplot(311)
p1, = ax2.plot(time, dFF, linewidth=2, color='green', label='dLight')
p2, = ax2.plot(LIGHT_x, y_scale*LIGHT_y+y_shift, linewidth=2, color='dodgerblue', label='Lick Event')
ax2.set_ylabel(r'$\Delta$F/F')
ax2.set_xlabel('Seconds')
ax2.set_title('dFF')
ax2.legend(handles=[p1], loc='upper right')
fig3.tight_layout()

#%% In [8] Modified z-score

#Calculate modified z-scores
z_mod_all = []
med_diff_all = []
mean_diff_all = []
med = np.median(Y_dF_all)
mean = np.mean(Y_dF_all)



for dF in dFF :
    med_diff = abs(dF-med)
    mean_diff = abs(dF-mean)
    med_diff_all.append(med_diff)
    mean_diff_all.append(mean_diff)


MedAD = np.median(med_diff_all)
MeanAD = np.mean(mean_diff_all)
mod_zmean = []
mod_zmed = []

for dF in dFF:
    mod_zmean.append((dF-med)/(1.253314*MeanAD))
    mod_zmed.append((dF-med)/(1.486*MedAD))

   
 
#%% Plot z-scored data
fig4 = plt.figure(figsize=(10,10))
ax3 = fig4.add_subplot(311)   

p2, = ax3.plot(time, mod_zmed, linewidth=2, color='green', label='dLight')
ax3.set_ylabel(r'Modified z-score')
ax3.set_xlabel('Seconds')
ax3.set_title('dFF z-scored')
ax3.legend(handles=[p1], loc='upper right')
fig4.tight_layout()

#%% Test of LIGHT_x (this doesn't really do anything)

# =============================================================================
# testarray =  [LIGHT_on, LIGHT_off]
# test2 = np.array([[1], [1]])
# test3 = np.kron(testarray, test2)
# test4 = test3.T
# test5 = np.reshape(test4, [1, -1])
# test6 = np.append(time[0], test5)
# test7 = np.append(time[0], test5[0])
# test8 = np.append(test7, time[-1])
# =============================================================================

#%% Average epoch data (wip - currently finding the time values for each individual epoch)

#pre_on = 1
#post_on = 1 #number of seconds after onset of stim delivery

#ttl_inds = []

#for i in LIGHT_on:
   # ttl_inds = np.where((time > LIGHT_on(i) - pre_on) & (time < LIGHT_on(i) + post_on))
 





