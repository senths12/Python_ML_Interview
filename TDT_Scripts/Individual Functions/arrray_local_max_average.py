import numpy as np
import matplotlib.pyplot as plt

time = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]) #time vector array
signal =  np.random.rand(15, 20000) #signal vector array
num_rows = np.shape(signal)[0] #number of rows for signal vector
Index_of_TTL = np.where(abs(time) == np.min(abs(time))) #for the 0, enter in the position of the x-coordinate of TTL signal
Index_of_TTL = int(Index_of_TTL[0]) #converts the 'tuple' of 1 element that is returned to an integer
Index_Final = np.where(abs(time-3)==np.min(abs(time-3))) #for the 5, enter in the position of the final time stamp
Index_Final = int(Index_Final[0]) #converts the 'tuple' of 1 element that is returned to an integer
rngArray = signal[0:num_rows, Index_of_TTL:Index_Final] #creates a new array (subarray) from the original array
localMax_percolumn = np.max(rngArray,axis=1) #find the maximum of each row in the subarray
trials = np.arange(1, (num_rows +1), 1) #creates the x-axis evenly spaced by 1 --> for the # of trials

ax1 = fig.add_subplot(312)
ax1.set_ylim(0, 4)
ax1.set_xlim(0, 15)
scatterplot = ax1.scatter(trials, localMax_percolumn) #creates a scatter plot between trials and 