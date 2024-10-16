import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def bin_array(arr, bin_size):
    # Convert to numpy array with float type for precision
    arr = np.array(arr, dtype=float)
    
    # Check if the array length is less than bin_size
    if len(arr) < bin_size:
        raise ValueError("Array length must be at least as large as the bin size.")
    
    # Calculate the number of bins
    num_bins = int(np.ceil(len(arr) / bin_size))
    
    # Pad the array to make sure it fits exactly into bins
    padded_length = num_bins * bin_size
    padded_arr = np.pad(arr, (0, padded_length - len(arr)), mode='constant', constant_values=np.nan)
    
    # Reshape and compute the mean of each bin
    reshaped_arr = padded_arr.reshape(num_bins, bin_size)
    binned_arr = np.nanmean(reshaped_arr, axis=1)
    
    return np.array(binned_arr)

def rms(arr, mean):
    return np.sqrt((np.nansum((arr - mean)**2))/len(arr))

path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
directory = os.fsencode(path_october)
#time of eclipse: [0:85], matching hours after eclipse: [66:245] (18:04-22:16)
RV_array = []
# nxt_day_time = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fits"):
        inputSpectrum = fits.open('{}/{}'.format(path_october, filename))
        RV_array.append(inputSpectrum[12].header["CCFRVMOD"] * 1000)
        # nxt_day_time.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"))

bin_arr = np.linspace(1,30,30)

# Bin the time series data
RMS_array = []
for i in range(0, len(bin_arr)):
    binned_data = bin_array(RV_array[66:245], int(bin_arr[i]))
    RMS_array.append(rms(binned_data, np.mean(binned_data)))

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
rv_obs = np.array(data["ccfrvmod"][130:-10]*1000 + 644.9)
# eclipse_time = []
# for i in data["obsdate"][130:-10]:
#     eclipse_time.append(datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5))

RMS_array_eclipse_day = []
for i in range(0, len(bin_arr)):
    binned_data = bin_array(rv_obs, int(bin_arr[i]))
    RMS_array_eclipse_day.append(rms(binned_data, np.mean(binned_data)))

# print(len(nxt_day_time[66:245]), len(eclipse_time))
# print(nxt_day_time[66], nxt_day_time[245], eclipse_time[0], eclipse_time[-1])
# 179 183
# 2023-10-15 18:04:59.501000 2023-10-15 22:16:15.170000 2023-10-14 18:05:13.500000 2023-10-14 22:17:28.500000

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(bin_arr, RMS_array, label = "10/15")
ax1.scatter(bin_arr, RMS_array_eclipse_day, label = "10/14 - eclipse")
ax1.set_xlabel("bin size")
ax1.set_ylabel("RMS (m/s)") 
plt.legend()
plt.savefig("photon_noise_binning.png")
plt.clf()