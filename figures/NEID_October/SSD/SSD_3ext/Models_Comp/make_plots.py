import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel

# GRASS
grass_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_all_lines_rv_regular_SSD_3ext.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
rv_error_GRASS_cb  = grass_data["rv_error"][()]
grass_data_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_all_lines_rv_off_SSD_3ext.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
rv_error_GRASS_no_cb = grass_data_no_cb["rv_error"][()]
# model
file = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_october_N_50_SSD_3ext.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]

# data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
rv_obs = np.array(rv_obs[0:-28])
rv_obs -= rv_obs[-1]
UTC_time = UTC_time[0:-28]

# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_lines = line_data["name"][()]
line_rv  = line_data["rv"][()]
rv_error_line = line_data["rv_error"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-28]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_line(UTC_time, line_rv_array_7, line_rv_array_8, label1, label2, save):
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, line_rv_array_7, color = 'k', marker = "x", s = 18, label = label1)
    axs[0].scatter(UTC_time, line_rv_array_8, color = 'y', marker = "x", s = 18, label = label2) 
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, line_rv_array_7 - line_rv_array_8, color = 'k', marker = "x", s = 3) 
    rms_model_no_cb = round(np.sqrt((np.nansum((line_rv_array_7 - line_rv_array_8)**2))/len(line_rv_array_7)),2)
    axs[1].text(UTC_time[-23], -5, "RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array_7[-4:-1] - line_rv_array_8[-4:-1])**2))/len(line_rv_array_7[-4:-1])),2)
    axs[1].text(UTC_time[-23], 10, "Out of Transit RMS {}".format(rms_grass_no_cb))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("{}.png".format(save))
    plt.clf()

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


line_rv_array_18 = jld2_read(line_data, line_rv, vb, 18)
line_rv_array_19 = jld2_read(line_data, line_rv, vb, 19)
plot_line(UTC_time, line_rv_array_18, line_rv_array_19, lines[18], lines[19], "unbinned")
plot_line(range(0, len(bin_array(line_rv_array_18,5))), bin_array(line_rv_array_18,5), bin_array(line_rv_array_19,5), lines[18], lines[19], "binned")

rv_error_line_array18 = line_data[rv_error_line[18]][()][0:-28]
rv_error_line_array19 = line_data[rv_error_line[19]][()][0:-28]
plot_line(UTC_time, rv_error_line_array18, rv_error_line_array19, lines[18], lines[19], "unbinned_err")
plot_line(range(0, len(bin_array(rv_error_line_array18,5))), bin_array(rv_error_line_array18,5), bin_array(rv_error_line_array19,5), lines[18], lines[19], "binned_err")