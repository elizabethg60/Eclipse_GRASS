import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel

grass_data_no_cb = h5py.File("SSD_3ext/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_ext_1.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv_lsf"][()]
brightness = grass_data_no_cb["brightness"][()]

# data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_rv  = line_data["rv"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs[0:-28])
rv_obs -= rv_obs[-1]

UTC_time = UTC_time[0:-28]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-28]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def jld2_read_fine(jld2_file, variable, vb, index):
    # original = jld2_file[variable[index]][()]
    # intensity = jld2_file[brightness[0]][()]
    # intensity = intensity / np.max(intensity)
    # array = [np.sum(original[i:i+56]*(intensity[i:i+56]))/np.sum(intensity[i:i+56]) for i in range(0, len(original), 56)]
    
    array = jld2_file[variable[index]][()][0:130]

    array = np.array(array + vb)
    array -= array[-1]
    return array

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

def bin_datetimes(arr, bin_size):
    # Convert datetime array to timestamps (float type) for precision
    arr = np.array([dt.timestamp() for dt in arr], dtype=float)
    
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
    
    # Convert back to datetime objects from the mean timestamp
    binned_datetimes = [datetime.fromtimestamp(ts) for ts in binned_arr]
    
    return np.array(binned_datetimes)

def plot_line(UTC_time, rv_obs, line_rv_array, GRASS, GRASS_label, path, save):
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'y', marker = "x", s = 18, label = "NEID RVs")
    axs[0].scatter(UTC_time, line_rv_array, color = 'k', marker = "x", s = 18, label = "NEID line RVs") 
    axs[0].plot(UTC_time, GRASS, color = 'b', linewidth = 2, label = GRASS_label)
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS)**2))/len(line_rv_array - GRASS)),2)
    print(rms_grass_no_cb)
    # axs[0].text(UTC_time[-40], 700, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[0].set_title("Fe I 5379 \u212B")
    # residuals
    axs[1].scatter(UTC_time, line_rv_array - GRASS, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("LunchTalk_KTS_figures/figure_6_{}.png".format(save))
    plt.clf()

for i in range(2,3):
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    GRASS_no_cb_array = jld2_read_fine(grass_data_no_cb, GRASS_no_cb, vb, i)

    # rm curve 
    # plot_line(UTC_time, rv_obs, line_rv_array, GRASS_no_cb_array, "GRASS Model", "No_Granulation/LinebyLine", lines_no_cb[i])
    plot_line(bin_datetimes(UTC_time,6), rv_obs, bin_array(line_rv_array,6), bin_array(GRASS_no_cb_array,6), "GRASS Model", "No_Granulation/LinebyLine", lines_no_cb[i])
