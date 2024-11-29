import os
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from astropy.time import Time
from barycorrpy import get_BC_vel
import matplotlib.colors as mcolors

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

# GRASS
grass_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_all_lines_rv_regular_SSD_3ext.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_rv  = line_data["rv"][()]
# data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[48:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

UTC_time = UTC_time[48:-28]
print(len(UTC_time))

path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
directory = os.fsencode(path_october)
RV_array_15 = []
nxt_day_time = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fits"):
        inputSpectrum = fits.open('{}/{}'.format(path_october, filename))
        RV_array_15.append(inputSpectrum[12].header["CCFRVMOD"] * 1000)
        nxt_day_time.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"))
# 0-83 roughly 2023-10-14 16:32:53.500000 2023-10-14 18:28:38.500000

path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/16/"
directory = os.fsencode(path_october)
RV_array_16 = []
nxt_day_time = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".fits"):
        inputSpectrum = fits.open('{}/{}'.format(path_october, filename))
        RV_array_16.append(inputSpectrum[12].header["CCFRVMOD"] * 1000)
        nxt_day_time.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"))

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][48:-28]
    array = np.array(array + vb)
    array -= array[-1]
    return array

bin_arr = np.linspace(1,15,15)

fig = plt.figure()
ax1 = fig.add_subplot()

# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=np.min(lines), vmax=np.max(lines) + 100)
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

for i, value in enumerate(lines):
    if i == 10:
        continue
    color = cmap(norm(value))

    GRASS_rv_array = jld2_read(grass_data, GRASS_rv, vb, i)
    line_rv_array = jld2_read(line_data, line_rv, vb, i)

    epsilon = 1e-10
    GRASS_rv_array_safe = np.where(GRASS_rv_array != 0, GRASS_rv_array, epsilon)
    residuals = line_rv_array / GRASS_rv_array_safe

    RMS_array_eclipse_day = []
    for j in range(0, len(bin_arr)):
        binned_data = bin_array(residuals, int(bin_arr[j]))
        RMS_array_eclipse_day.append(rms(binned_data, np.mean(binned_data)))
    ax1.scatter(bin_arr, RMS_array_eclipse_day, color = color)
    ax1.plot(bin_arr, RMS_array_eclipse_day, color = color)

RMS_array_15 = []
for i in range(0, len(bin_arr)):
        binned_data = bin_array(RV_array_15[0:83], int(bin_arr[i]))
        RMS_array_15.append(rms(binned_data, np.mean(binned_data)))
ax1.scatter(bin_arr, RMS_array_15, color = 'g', label = "10/15")
ax1.plot(bin_arr, RMS_array_15, color = 'g')

RMS_array_16 = []
for i in range(0, len(bin_arr)):
        binned_data = bin_array(RV_array_16[0:83], int(bin_arr[i]))
        RMS_array_16.append(rms(binned_data, np.mean(binned_data)))
ax1.scatter(bin_arr, RMS_array_16, color = 'y', label = "10/16")
ax1.plot(bin_arr, RMS_array_16, color = 'y')

ax1.set_xlabel("bin size")
ax1.set_ylabel("RMS (m/s)") 
plt.legend()
plt.yscale('log')
plt.savefig("noise_binning.png")
plt.clf()

