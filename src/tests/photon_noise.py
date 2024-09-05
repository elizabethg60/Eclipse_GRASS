import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# def bin_array(arr, bin_size):
#     # Convert to numpy array with float type for precision
#     arr = np.array(arr, dtype=float)
    
#     # Check if the array length is less than bin_size
#     if len(arr) < bin_size:
#         raise ValueError("Array length must be at least as large as the bin size.")
    
#     # Calculate the number of bins
#     num_bins = int(np.ceil(len(arr) / bin_size))
    
#     # Pad the array to make sure it fits exactly into bins
#     padded_length = num_bins * bin_size
#     padded_arr = np.pad(arr, (0, padded_length - len(arr)), mode='constant', constant_values=np.nan)
    
#     # Reshape and compute the mean of each bin
#     reshaped_arr = padded_arr.reshape(num_bins, bin_size)
#     binned_arr = np.nanmean(reshaped_arr, axis=1)
    
#     return np.array(binned_arr)

import pandas as pd
import numpy as np

def bin_time_series(timestamps, velocities, bin_intervals):
    # Convert timestamps to pandas datetime format
    timestamps = pd.to_datetime(timestamps)
    
    # Create a DataFrame from the velocities and timestamps
    df = pd.DataFrame({
        'timestamp': timestamps,
        'velocity': velocities
    })
    
    # Set the timestamp column as the index
    df.set_index('timestamp', inplace=True)
    
    # Dictionary to hold binned data for each interval
    binned_data = {}
    
    for interval in bin_intervals:
        # Resample the data to the specified interval
        resampled_df = df.resample(interval).agg({'velocity': ['mean', 'std']})
        
        # Rename columns for clarity
        resampled_df.columns = ['mean_velocity', 'velocity_error']
        
        # Store the result in the dictionary
        binned_data[interval] = resampled_df
        
    return binned_data

def rms(arr, mean):
    return np.sqrt((np.nansum((arr - mean)**2))/len(arr))

path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
directory = os.fsencode(path_october)

#time of eclipse: [0:85], matching hours after eclipse: [66:245] (18:04-22:16)
RV_array = []
nxt_day_time = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fits"):
        inputSpectrum = fits.open('{}/{}'.format(path_october, filename))
        RV_array.append(inputSpectrum[12].header["CCFRVMOD"] * 1000)
        nxt_day_time.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"))

# RMS_array = []
    # Define the bin intervals you want (e.g., 1min, 5min, 10min)
bin_intervals = ['1T', '2T', '3T', '4T', '5T', '6T', '7T', '8T', '9T', '10T', '11T', '12T', '13T','14T', '15T', '16T', '17T', '18T','19T', '20T', '21T', '22T', '23T', '24T', '25T', '26T', '27T', '28T', '29T', '30T']
    
# Bin the time series data
RMS_array = []
err_nxt_day = []
bin_arr = np.linspace(1,30,30)
for i in range(0, len(bin_intervals)):
    binned_data = bin_time_series(nxt_day_time[66:245], RV_array[66:245], [bin_intervals[i]])
    RMS_array.append(rms(binned_data[bin_intervals[i]]['mean_velocity'], np.mean(RV_array[66:245])))
    err_nxt_day.append(np.mean(binned_data[bin_intervals[i]]['velocity_error']))

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
rv_obs = np.array(data["ccfrvmod"][130:-10]*1000 + 644.9)
eclipse_time = []
for i in data["obsdate"][130:-10]:
    eclipse_time.append(datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5))

RMS_array_eclipse_day = []
err_eclipse = []
for i in range(0, len(bin_intervals)):
    binned_data = bin_time_series(eclipse_time, rv_obs, [bin_intervals[i]])
    RMS_array_eclipse_day.append(rms(binned_data[bin_intervals[i]]['mean_velocity'], np.mean(rv_obs)))
    err_eclipse.append(np.mean(binned_data[bin_intervals[i]]['velocity_error']))

# print(len(nxt_day_time[66:245]), len(eclipse_time))
# print(nxt_day_time[66], nxt_day_time[245], eclipse_time[0], eclipse_time[-1])
# 179 183
# 2023-10-15 18:04:59.501000 2023-10-15 22:16:15.170000 2023-10-14 18:05:13.500000 2023-10-14 22:17:28.500000

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.errorbar(bin_arr, RMS_array, label = "10/15", yerr=err_nxt_day, fmt='o')
ax1.errorbar(bin_arr, RMS_array_eclipse_day, label = "10/14 - eclipse", yerr=err_eclipse, fmt='o')
ax1.set_xlabel("bin size (minutes)")
ax1.set_ylabel("RMS (m/s)") 
plt.legend()
plt.savefig("photon_noise_binning.png")
plt.clf()