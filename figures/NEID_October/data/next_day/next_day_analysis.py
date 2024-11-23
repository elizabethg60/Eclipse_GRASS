import h5py
import numpy as np
import pandas as pd
import csv
from astropy.time import Time
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

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

def rms(arr):
    return np.sqrt((np.nansum((arr - np.mean(arr))**2))/len(arr)) 

neid_1015 = ["2023-10-15T16:33:02.500", "2023-10-15T16:34:25.500", "2023-10-15T16:35:48.500", "2023-10-15T16:37:10.500", "2023-10-15T16:38:33.500", "2023-10-15T16:39:56.500", "2023-10-15T16:41:18.500", "2023-10-15T16:42:41.500", "2023-10-15T16:44:04.500", "2023-10-15T16:45:26.500", "2023-10-15T16:46:49.500", "2023-10-15T16:48:12.500", "2023-10-15T16:49:34.500", "2023-10-15T16:50:57.500", "2023-10-15T16:52:20.500", "2023-10-15T16:53:42.500", "2023-10-15T16:55:05.500", "2023-10-15T16:56:28.500", "2023-10-15T16:57:50.500", "2023-10-15T16:59:13.500", "2023-10-15T17:00:36.500", "2023-10-15T17:01:59.500", "2023-10-15T17:03:21.500", "2023-10-15T17:04:44.500", "2023-10-15T17:06:07.500", "2023-10-15T17:07:29.500", "2023-10-15T17:08:52.500", "2023-10-15T17:10:15.500", "2023-10-15T17:11:37.500", "2023-10-15T17:13:00.500", "2023-10-15T17:14:23.500", "2023-10-15T17:15:45.500", "2023-10-15T17:17:08.500", "2023-10-15T17:18:31.500", "2023-10-15T17:19:53.500", "2023-10-15T17:21:16.500", "2023-10-15T17:22:39.500", "2023-10-15T17:24:01.500", "2023-10-15T17:25:24.500", "2023-10-15T17:26:47.500", "2023-10-15T17:28:09.500", "2023-10-15T17:29:32.500", "2023-10-15T17:30:55.500", "2023-10-15T17:32:18.500", "2023-10-15T17:33:40.500", "2023-10-15T17:35:03.500", "2023-10-15T17:36:26.500", "2023-10-15T17:37:48.500", "2023-10-15T17:39:11.500", "2023-10-15T17:40:34.500", "2023-10-15T17:41:56.500", "2023-10-15T17:43:19.500", "2023-10-15T17:44:42.500", "2023-10-15T17:46:04.500", "2023-10-15T17:47:27.500", "2023-10-15T17:48:50.500", "2023-10-15T17:50:12.500", "2023-10-15T17:51:35.500", "2023-10-15T17:52:58.500", "2023-10-15T17:54:20.500", "2023-10-15T17:57:10.500", "2023-10-15T17:58:33.500", "2023-10-15T17:59:55.500", "2023-10-15T18:01:18.500", "2023-10-15T18:02:41.500", "2023-10-15T18:04:03.500", "2023-10-15T18:05:26.500", "2023-10-15T18:06:49.500", "2023-10-15T18:08:11.500", "2023-10-15T18:09:34.500", "2023-10-15T18:10:57.500", "2023-10-15T18:12:19.500", "2023-10-15T18:13:42.500", "2023-10-15T18:15:05.500", "2023-10-15T18:16:27.500", "2023-10-15T18:17:50.500", "2023-10-15T18:19:13.500", "2023-10-15T18:20:35.500", "2023-10-15T18:21:58.500", "2023-10-15T18:23:21.500", "2023-10-15T18:24:44.500", "2023-10-15T18:26:06.500", "2023-10-15T18:27:29.500", "2023-10-15T18:28:52.500", "2023-10-15T18:30:14.500", "2023-10-15T18:31:37.500", "2023-10-15T18:33:00.500", "2023-10-15T18:34:22.500", "2023-10-15T18:35:45.500", "2023-10-15T18:37:08.500", "2023-10-15T18:38:30.500", "2023-10-15T18:39:53.500", "2023-10-15T18:41:16.500", "2023-10-15T18:42:38.500", "2023-10-15T18:44:01.500", "2023-10-15T18:45:24.500", "2023-10-15T18:46:46.500", "2023-10-15T18:48:09.500", "2023-10-15T18:49:32.500", "2023-10-15T18:50:54.500", "2023-10-15T18:52:17.500", "2023-10-15T18:53:40.500", "2023-10-15T18:55:03.500", "2023-10-15T18:56:25.500", "2023-10-15T18:57:48.500", "2023-10-15T18:59:11.500", "2023-10-15T19:00:33.500", "2023-10-15T19:01:56.500", "2023-10-15T19:03:19.500", "2023-10-15T19:04:41.500", "2023-10-15T19:06:04.500"]
time_julian = []
for i in neid_1015:
    dt = datetime.strptime(i, "%Y-%m-%dT%H:%M:%S.%f")
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

#line by line data
line_data = h5py.File("neid_RVlinebyline_nxt_day.jld2", "r")
lines = line_data["name"][()]
line_rv = line_data["rv"][()]
line_rv_err = line_data["rv_error"][()]

#projected rv model - regular
file_regular_SSD = h5py.File("neid_october_N_50_nxt_day_SSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD = h5py.File("neid_all_lines_rv_regular_SSD_nxt_day.jld2", "r")
GRASS_rv_KSSD  = grass_data_KSSD["rv"][()]
rv_error_GRASS_cb  = grass_data_KSSD["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD = h5py.File("neid_all_lines_rv_off_SSD_nxt_day.jld2", "r")
GRASS_no_cb_v_KSSD  = grass_data_no_cb_KSSD["rv"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_3(projected, grass_cb, grass_no_cb, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, title, line_file, line_rv,
                 label1, label2, label3):
    
    projected_arr_pipe = []
    line_rv_pipe = []
    grass_cb_arr_pipe = []
    grass_no_cb_arr_pipe = []
    lines_arr = []
    line_err = []
    lines_combined = []
    for i in range(0,len(lines)):
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-25]))

        lines_arr.append(lines[i])

        projected_arr_pipe.append(rms(bin_array(projected_i,6)))
        grass_cb_arr_pipe.append(rms(bin_array(grass_cb_i,6)))
        grass_no_cb_arr_pipe.append(rms(bin_array(grass_no_cb_i,6)))
        line_rv_pipe.append(rms(bin_array(line_i,6)))

        lines_combined.append(list(line_i-grass_no_cb_i))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none', label = label1)
    plt.scatter(lines_arr, grass_cb_arr_pipe, color = 'r', facecolors='none', label = label2)
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none', label = label3)
    plt.scatter(lines_arr, line_rv_pipe, color = 'y', facecolors='none', label = "line rv")
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

    lines_combined_avg = np.mean(lines_combined, axis=0)
    print(rms(bin_array(lines_combined_avg,6)))

plot_3(RV_list_no_cb_SSD, GRASS_rv_KSSD, GRASS_no_cb_v_KSSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_KSSD, grass_data_no_cb_KSSD, "next_day_analysis", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
