import h5py
import numpy as np
import pandas as pd
import csv
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from barycorrpy import get_BC_vel

def bin_array(arr, bin_size):
# bin_size = 5 ~ 6 minute bins
# cloudless data = 130 points so when binned by 5 get 26 binned points each with 5 points
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

def rms(arr, model):
    return np.sqrt((np.nansum((arr - model)**2))/len(arr))   

def rms_single(arr):
    return np.sqrt((np.nansum((arr)**2))/len(arr))    

def weighted_rms(arr, model, weight):
    return np.sqrt((np.nansum((arr - model)**2 * weight))/(np.nansum(weight)*len(arr)))

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
time_julian = []
UTC_time = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)[0:-28]
UTC_time = UTC_time[0:-28]

# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
lines = line_data["name"][()]
line_rv = line_data["rv"][()]
line_rv_err = line_data["rv_error"][()]

# NL94
# projected rv model - regular
file_regular_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_october_N_50_NL94.jld2", "r")
RV_list_no_cb_NL94 = file_regular_NL94["RV_list_no_cb"][()]
# GRASS CB
grass_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_all_lines_rv_regular_NL94.jld2", "r")
GRASS_rv_NL94  = grass_data_NL94["rv"][()]
# GRASS no CB
grass_data_no_cb_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_all_lines_rv_off_NL94.jld2", "r")
GRASS_no_cb_v_NL94  = grass_data_no_cb_NL94["rv"][()]

# 300
file_regular_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300/data/neid_october_N_50_300.jld2", "r")
RV_list_no_cb_300 = file_regular_300["RV_list_no_cb"][()]
# GRASS CB
grass_data_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300/data/neid_all_lines_rv_regular_300.jld2", "r")
GRASS_rv_300  = grass_data_300["rv"][()]
# GRASS no CB
grass_data_no_cb_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300/data/neid_all_lines_rv_off_300.jld2", "r")
GRASS_no_cb_v_300  = grass_data_no_cb_300["rv"][()]

# 300 + 3ext
file_regular_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300_3ext/data/neid_october_N_50_300_3ext.jld2", "r")
RV_list_no_cb_300_3ext = file_regular_300_3ext["RV_list_no_cb"][()]
# GRASS CB
grass_data_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300_3ext/data/neid_all_lines_rv_regular_300_3ext.jld2", "r")
GRASS_rv_300_3ext  = grass_data_300_3ext["rv"][()]
# GRASS no CB
grass_data_no_cb_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/300/300_3ext/data/neid_all_lines_rv_off_300_3ext.jld2", "r")
GRASS_no_cb_v_300_3ext  = grass_data_no_cb_300_3ext["rv"][()]

# SSD
# projected rv model - regular
file_regular_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD/data/neid_october_N_50_SSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
# GRASS CB
grass_data_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD/data/neid_all_lines_rv_regular_SSD.jld2", "r")
GRASS_rv_SSD  = grass_data_SSD["rv"][()]
rv_error_GRASS_cb  = grass_data_SSD["rv_error"][()]
# GRASS no CB
grass_data_no_cb_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD/data/neid_all_lines_rv_off_SSD.jld2", "r")
GRASS_no_cb_v_SSD  = grass_data_no_cb_SSD["rv"][()]

# SSD + 3ext
# projected rv model - regular
file_regular_SSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_october_N_50_SSD_3ext.jld2", "r")
RV_list_no_cb_SSD_3ext = file_regular_SSD_3ext["RV_list_no_cb"][()]
# GRASS CB
grass_data_SSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_all_lines_rv_regular_SSD_3ext.jld2", "r")
GRASS_rv_SSD_3ext  = grass_data_SSD_3ext["rv"][()]
rv_error_GRASS_cb_3ext  = grass_data_SSD_3ext["rv_error"][()]
# GRASS no CB
grass_data_no_cb_SSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD_3ext/data/neid_all_lines_rv_off_SSD_3ext.jld2", "r")
GRASS_no_cb_v_SSD_3ext  = grass_data_no_cb_SSD_3ext["rv"][()]

# HD
# projected rv model - regular
file_regular_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/HD/data/neid_october_N_50_HD.jld2", "r")
RV_list_no_cb_HD = file_regular_HD["RV_list_no_cb"][()]
# GRASS CB
grass_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/HD/data/neid_all_lines_rv_regular_HD.jld2", "r")
GRASS_rv_HD  = grass_data_HD["rv"][()]
# GRASS no CB
grass_data_no_cb_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/HD/data/neid_all_lines_rv_off_HD.jld2", "r")
GRASS_no_cb_v_HD  = grass_data_no_cb_HD["rv"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-28]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_line_bin(projected, grass_cb, grass_no_cb, 
                 projected_file, grass_cb_file, grass_no_cb_file, title, line_file, line_rv,
                 label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    lines_arr = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(bin_array(line_i,5), bin_array(projected_i,5)))
        grass_cb_arr.append(rms(bin_array(line_i,5), bin_array(grass_cb_i,5)))
        grass_no_cb_arr.append(rms(bin_array(line_i,5), bin_array(grass_no_cb_i,5)))

        # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
        # axs[0].set_title("Projected RV (no granulation mitigation)")
        # axs[0].scatter(UTC_time, line_i, color = 'k', marker = "x", s = 18)
        # axs[0].plot(UTC_time, projected_i, color = 'r', linewidth = 2)
        # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[0].set_xlabel("Time (UTC)", fontsize=12)
        # axs[0].set_ylabel("RV [m/s]", fontsize=12)
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # # residuals
        # axs[1].scatter(UTC_time, line_i - projected_i, color = 'r', marker = "x", s = 3) 
        # rms_model_no_cb = rms(bin_array(line_i,5), bin_array(projected_i,5))
        # axs[1].text(UTC_time[-30], 30, "RMS {} m/s".format(round(rms_model_no_cb,2)))
        # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[1].set_xlabel("Time (UTC)", fontsize=12)
        # axs[1].set_ylabel("Residuals", fontsize=12) 
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # plt.savefig("model1.pdf")
        # plt.clf()

        # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
        # axs[0].set_title("GRASS (no temporal variability)")
        # axs[0].scatter(UTC_time, line_i, color = 'k', marker = "x", s = 18)
        # axs[0].plot(UTC_time, grass_no_cb_i, color = 'r', linewidth = 2)
        # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[0].set_xlabel("Time (UTC)", fontsize=12)
        # axs[0].set_ylabel("RV [m/s]", fontsize=12)
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # # residuals
        # axs[1].scatter(UTC_time, line_i - grass_no_cb_i, color = 'r', marker = "x", s = 3) 
        # rms_model_no_cb = rms(bin_array(line_i,5), bin_array(grass_no_cb_i,5))
        # axs[1].text(UTC_time[-30], 30, "RMS {} m/s".format(round(rms_model_no_cb,2)))
        # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[1].set_xlabel("Time (UTC)", fontsize=12)
        # axs[1].set_ylabel("Residuals", fontsize=12) 
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # plt.savefig("model2.pdf")
        # plt.clf()

        # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
        # axs[0].set_title("GRASS")
        # axs[0].scatter(UTC_time, line_i, color = 'k', marker = "x", s = 18)
        # axs[0].plot(UTC_time, grass_cb_i, color = 'r', linewidth = 2)
        # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[0].set_xlabel("Time (UTC)", fontsize=12)
        # axs[0].set_ylabel("RV [m/s]", fontsize=12)
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # # residuals
        # axs[1].scatter(UTC_time, line_i - grass_cb_i, color = 'r', marker = "x", s = 3) 
        # rms_model_no_cb = rms(bin_array(line_i,5), bin_array(grass_cb_i,5))
        # axs[1].text(UTC_time[-30], 30, "RMS {} m/s".format(round(rms_model_no_cb,2)))
        # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        # axs[1].set_xlabel("Time (UTC)", fontsize=12)
        # axs[1].set_ylabel("Residuals", fontsize=12) 
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)
        # plt.savefig("model3.pdf")
        # plt.clf()

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_model_bin(projected, grass_cb, grass_no_cb, data_300, 
                 projected_file, grass_cb_file, grass_no_cb_file, file_300, title, 
                 line_file, line_rv, label1, label2, label3, label4):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    lines_arr = []
    arr_300 = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        data_300_i = jld2_read(file_300, data_300, vb, i)
        
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(bin_array(line_i,5), bin_array(projected_i,5)))
        grass_cb_arr.append(rms(bin_array(line_i,5), bin_array(grass_cb_i,5)))
        grass_no_cb_arr.append(rms(bin_array(line_i,5), bin_array(grass_no_cb_i,5)))
        arr_300.append(rms(bin_array(line_i,5), bin_array(data_300_i,5)))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, arr_300, label = label4, color = 'y')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_3(projected, grass_cb, grass_no_cb, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, title, line_file, line_rv,
                 label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    projected_arr_pipe = []
    grass_cb_arr_pipe = []
    grass_no_cb_arr_pipe = []
    projected_arr_weight = []
    grass_cb_arr_weight = []
    grass_no_cb_arr_weight = []
    lines_arr = []
    grass_cb_err_arr = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))

        projected_arr_weight.append(weighted_rms(line_i, projected_i, line_file[line_rv_err[i]][()][0:-28]))
        grass_cb_arr_weight.append(weighted_rms(line_i, grass_cb_i, line_file[line_rv_err[i]][()][0:-28]))
        grass_no_cb_arr_weight.append(weighted_rms(line_i, grass_no_cb_i, line_file[line_rv_err[i]][()][0:-28]))

        if title == "SSD" or title == "LD_comp":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-28]))
    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_pipe, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none')
    plt.scatter(lines_arr, projected_arr_weight, color = 'b', marker='x')
    plt.scatter(lines_arr, grass_cb_arr_weight, color = 'r', marker='x')
    plt.scatter(lines_arr, grass_no_cb_arr_weight, color = 'g', marker='x')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')
    if title == "SSD" or title == "LD_comp":
        plt.scatter(lines_arr, grass_cb_err_arr, color = 'k', marker = 'x')
        plt.title("solid - line rms & circle - pipeline rms & x - mean error grass cb ccf (SSD)")
    else:
        plt.title("solid - line rms & circle - pipeline rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_4_4(projected, grass_cb, grass_no_cb, data_300, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, file_300, title, 
                 line_file, line_rv, label1, label2, label3, label4):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    projected_arr_pipe = []
    grass_cb_arr_pipe = []
    grass_no_cb_arr_pipe = []
    lines_arr = []
    grass_cb_err_arr = []
    arr_300 = []
    arr_300_pipe = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        data_300_i = jld2_read(file_300, data_300, vb, i)
        
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))
        arr_300.append(rms(line_i, data_300_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))
        arr_300_pipe.append(rms(rv_obs, data_300_i))

        if title == "SSD" or title == "LD_comp":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-28]))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, arr_300, label = label4, color = 'k')
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_pipe, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none')
    plt.scatter(lines_arr, arr_300_pipe, color = 'k', facecolors='none')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')

    if title == "SSD" or title == "LD_comp":
        plt.scatter(lines_arr, grass_cb_err_arr, color = 'k', marker = 'x')
        plt.title("solid - line rms & circle - pipeline rms & x - mean error grass cb ccf (SSD)")
    else:
        plt.title("solid - line rms & circle - pipeline rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_2_2(projected, grass_no_cb, grass_cb_err, 
                 projected_file, grass_no_cb_file, title, 
                 line_file, line_rv, label1, label3):
    
    projected_arr = []
    grass_no_cb_arr = []
    projected_arr_pipe = []
    grass_no_cb_arr_pipe = []
    lines_arr = []
    grass_cb_err_arr = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')
    plt.title("solid - line rms & circle - pipeline rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf()  

def plot_bin(projected, grass_no_cb, grass_cb_err, 
                 projected_file, grass_no_cb_file, title, line_file, line_rv,
                 label1, label3):
    
    projected_arr = []
    grass_no_cb_arr = []
    projected_arr_bin = []
    grass_no_cb_arr_bin = []
    lines_arr = []
    line_err = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_bin.append(rms(bin_array(line_i,5), bin_array(projected_i,5)))
        grass_no_cb_arr_bin.append(rms(bin_array(line_i,5), bin_array(grass_no_cb_i,5)))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_bin, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_bin, color = 'g', facecolors='none')
    plt.scatter(lines_arr, line_err, color = 'k', marker='x')

    plt.title("solid - line rms & circle - bin rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf()     

def plot_single_rms(line_file, line_rv, projected_file, projected, grass_cb_file, grass_cb, grass_no_cb_file, grass_no_cb, 
                    projected_file_3ext, projected_3ext, grass_cb_file_3ext, grass_cb_3ext, grass_no_cb_file_3ext, grass_no_cb_3ext,
                    title):
    SSD_projected = []
    SSD_grass_cb = []
    SSD_grass_no_cb = []
    SSD_3ext_projected = []
    SSD_3ext_grass_cb = []
    SSD_3ext_grass_no_cb = []
    SSD_projected_pipeline = []
    SSD_grass_cb_pipeline = []
    SSD_grass_no_cb_pipeline = []
    SSD_3ext_projected_pipeline = []
    SSD_3ext_grass_cb_pipeline = []
    SSD_3ext_grass_no_cb_pipeline = []

    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        projected_i_3ext = jld2_read(projected_file_3ext, projected_3ext, vb, i)
        grass_cb_i_3ext = jld2_read(grass_cb_file_3ext, grass_cb_3ext, vb, i)
        grass_no_cb_i_3ext = jld2_read(grass_no_cb_file_3ext, grass_no_cb_3ext, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)

        for j in range(0, len(line_i)):
            SSD_projected.append(line_i[j] - projected_i[j])
            SSD_grass_cb.append(line_i[j] - grass_cb_i[j])
            SSD_grass_no_cb.append(line_i[j] - grass_no_cb_i[j])
            SSD_3ext_projected.append(line_i[j] - projected_i_3ext[j])
            SSD_3ext_grass_cb.append(line_i[j] - grass_cb_i_3ext[j])
            SSD_3ext_grass_no_cb.append(line_i[j] - grass_no_cb_i_3ext[j])

            SSD_projected_pipeline.append(rv_obs[j] - projected_i[j])
            SSD_grass_cb_pipeline.append(rv_obs[j] - grass_cb_i[j])
            SSD_grass_no_cb_pipeline.append(rv_obs[j] - grass_no_cb_i[j])
            SSD_3ext_projected_pipeline.append(rv_obs[j] - projected_i_3ext[j])
            SSD_3ext_grass_cb_pipeline.append(rv_obs[j] - grass_cb_i_3ext[j])
            SSD_3ext_grass_no_cb_pipeline.append(rv_obs[j] - grass_no_cb_i_3ext[j])

    SSD_projected = np.array(SSD_projected)
    SSD_grass_cb = np.array(SSD_grass_cb)
    SSD_grass_no_cb = np.array(SSD_grass_no_cb)
    SSD_3ext_projected = np.array(SSD_3ext_projected)
    SSD_3ext_grass_cb = np.array(SSD_3ext_grass_cb)
    SSD_3ext_grass_no_cb = np.array(SSD_3ext_grass_no_cb)

    SSD_projected_pipeline = np.array(SSD_projected_pipeline)
    SSD_grass_cb_pipeline = np.array(SSD_grass_cb_pipeline)
    SSD_grass_no_cb_pipeline = np.array(SSD_grass_no_cb_pipeline)
    SSD_3ext_projected_pipeline = np.array(SSD_3ext_projected_pipeline)
    SSD_3ext_grass_cb_pipeline = np.array(SSD_3ext_grass_cb_pipeline)
    SSD_3ext_grass_no_cb_pipeline = np.array(SSD_3ext_grass_no_cb_pipeline)

    x = [0,1,2]   
    y_SSD = [np.sqrt((np.nansum((SSD_projected)**2))/len(SSD_projected)), np.sqrt((np.nansum((SSD_grass_cb)**2))/len(SSD_grass_cb)), np.sqrt((np.nansum((SSD_grass_no_cb)**2))/len(SSD_grass_no_cb))]
    y_SSD_3ext = [np.sqrt((np.nansum((SSD_3ext_projected)**2))/len(SSD_3ext_projected)), np.sqrt((np.nansum((SSD_3ext_grass_cb)**2))/len(SSD_3ext_grass_cb)), np.sqrt((np.nansum((SSD_3ext_grass_no_cb)**2))/len(SSD_3ext_grass_no_cb))]
    y_SSD_pipeline = [np.sqrt((np.nansum((SSD_projected_pipeline)**2))/len(SSD_projected_pipeline)), np.sqrt((np.nansum((SSD_grass_cb_pipeline)**2))/len(SSD_grass_cb_pipeline)), np.sqrt((np.nansum((SSD_grass_no_cb_pipeline)**2))/len(SSD_grass_no_cb_pipeline))]
    y_SSD_3ext_pipeline = [np.sqrt((np.nansum((SSD_3ext_projected_pipeline)**2))/len(SSD_3ext_projected_pipeline)), np.sqrt((np.nansum((SSD_3ext_grass_cb_pipeline)**2))/len(SSD_3ext_grass_cb_pipeline)), np.sqrt((np.nansum((SSD_3ext_grass_no_cb_pipeline)**2))/len(SSD_3ext_grass_no_cb_pipeline))]
    labels = ['Projected', 'GRASS Var', 'GRASS No Var']
    plt.figure(figsize=(12, 6))
    plt.scatter(x, y_SSD, color = 'b', label = "SSD")
    plt.scatter(x, y_SSD_3ext, color = 'g', label = "SSD+3ext")

    plt.scatter(x, y_SSD_pipeline,  color = 'b', facecolors='none')
    plt.scatter(x, y_SSD_3ext_pipeline, color = 'g', facecolors='none')
    plt.title("solid - line rms & circle - pipeline rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(x, labels)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.legend()
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def rms_comp(line_file, line_rv, title):
    line_err = []
    new_err = []
    for i in range(0,len(lines)):

        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err_i = line_file[line_rv_err[i]][()][0:-28]
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-28]))

        new_rms = []
        for j in range(0, len(line_err_i)):
            new_rms.append(np.random.normal(0.0, line_err_i[j]))
        new_err.append(rms_single(np.array(new_rms)))

    plt.scatter(lines, line_err,  color = 'b', label='Mean St. Error')
    plt.scatter(lines, new_err, color = 'g', label='Random Draw')
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.legend()
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

# plot_line_bin(RV_list_no_cb_NL94, GRASS_rv_NL94, GRASS_no_cb_v_NL94, file_regular_NL94, grass_data_NL94, grass_data_no_cb_NL94, "Bins/NL94", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_line_bin(RV_list_no_cb_SSD, GRASS_rv_SSD, GRASS_no_cb_v_SSD, file_regular_SSD, grass_data_SSD, grass_data_no_cb_SSD, "Bins/SSD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_line_bin(RV_list_no_cb_HD, GRASS_rv_HD, GRASS_no_cb_v_HD, file_regular_HD, grass_data_HD, grass_data_no_cb_HD, "Bins/HD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_line_bin(RV_list_no_cb_300, GRASS_rv_300, GRASS_no_cb_v_300, file_regular_300, grass_data_300, grass_data_no_cb_300, "Bins/300", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_line_bin(RV_list_no_cb_SSD_3ext, GRASS_rv_SSD_3ext, GRASS_no_cb_v_SSD_3ext, file_regular_SSD_3ext, grass_data_SSD_3ext, grass_data_no_cb_SSD_3ext, "Bins/SSD_3ext", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')

# plot_model_bin(RV_list_no_cb_NL94, RV_list_no_cb_SSD, RV_list_no_cb_HD, RV_list_no_cb_300, file_regular_NL94, file_regular_SSD, file_regular_HD, file_regular_300, "Bins/Projected", line_data, line_rv, "NL94", "SSD", "HD", "300")
# plot_model_bin(GRASS_rv_NL94, GRASS_rv_SSD, GRASS_rv_HD, GRASS_rv_300, grass_data_NL94, grass_data_SSD, grass_data_HD, grass_data_300, "Bins/GRASS_cb", line_data, line_rv, "NL94", "SSD", "HD", "300")
# plot_model_bin(GRASS_no_cb_v_NL94, GRASS_no_cb_v_SSD, GRASS_no_cb_v_HD, GRASS_no_cb_v_300, grass_data_no_cb_NL94, grass_data_no_cb_SSD, grass_data_no_cb_HD, grass_data_no_cb_300, "Bins/GRASS_no_cb", line_data, line_rv, "NL94", "SSD", "HD", "300")

# plot_3(RV_list_no_cb_NL94, GRASS_rv_NL94, GRASS_no_cb_v_NL94, 0, file_regular_NL94, grass_data_NL94, grass_data_no_cb_NL94, "LD_comp/NL94", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_3(RV_list_no_cb_SSD, GRASS_rv_SSD, GRASS_no_cb_v_SSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_SSD, grass_data_no_cb_SSD, "LD_comp/SSD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_3(RV_list_no_cb_HD, GRASS_rv_HD, GRASS_no_cb_v_HD, 0, file_regular_HD, grass_data_HD, grass_data_no_cb_HD, "LD_comp/HD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_3(RV_list_no_cb_300, GRASS_rv_300, GRASS_no_cb_v_300, 0, file_regular_300, grass_data_300, grass_data_no_cb_300, "LD_comp/300", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_4_4(GRASS_rv_NL94, GRASS_rv_SSD, GRASS_rv_HD, GRASS_rv_300, rv_error_GRASS_cb, grass_data_NL94, grass_data_SSD, grass_data_HD, grass_data_300, "LD_comp/LD_comp", line_data, line_rv, "NL94", "SSD", "HD", "300")

# plot_3(RV_list_no_cb_SSD, GRASS_rv_SSD, GRASS_no_cb_v_SSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_SSD, grass_data_no_cb_SSD, "SSD_comp/SSD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_3(RV_list_no_cb_SSD_3ext, GRASS_rv_SSD_3ext, GRASS_no_cb_v_SSD_3ext, rv_error_GRASS_cb_3ext, file_regular_SSD_3ext, grass_data_SSD_3ext, grass_data_no_cb_SSD_3ext, "SSD_comp/SSD_3ext", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_2_2(GRASS_rv_SSD, GRASS_rv_SSD_3ext, rv_error_GRASS_cb, grass_data_SSD, grass_data_SSD_3ext, "SSD_comp/grass_cb_vs_ext", line_data, line_rv, "SSD", "SSD+3ext")
# plot_bin(GRASS_rv_SSD, GRASS_rv_SSD_3ext, rv_error_GRASS_cb, grass_data_SSD, grass_data_SSD_3ext, "SSD_comp/grass_cb_vs_ext_binned", line_data, line_rv, "SSD", "SSD+3ext")
# plot_single_rms(line_data, line_rv, file_regular_SSD, RV_list_no_cb_SSD, grass_data_SSD, GRASS_rv_SSD, grass_data_no_cb_SSD, GRASS_no_cb_v_SSD,
#                 file_regular_SSD_3ext, RV_list_no_cb_SSD_3ext, grass_data_SSD_3ext, GRASS_rv_SSD_3ext, grass_data_no_cb_SSD_3ext, GRASS_no_cb_v_SSD_3ext,
#                 "SSD_comp/single_rms")

# plot_3(RV_list_no_cb_300, GRASS_rv_300, GRASS_no_cb_v_300, rv_error_GRASS_cb, file_regular_300, grass_data_300, grass_data_no_cb_300, "300_comp/300", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_3(RV_list_no_cb_300_3ext, GRASS_rv_300_3ext, GRASS_no_cb_v_300_3ext, rv_error_GRASS_cb_3ext, file_regular_300_3ext, grass_data_300_3ext, grass_data_no_cb_300_3ext, "300_comp/300_3ext", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
# plot_2_2(GRASS_rv_300, GRASS_rv_300_3ext, rv_error_GRASS_cb, grass_data_300, grass_data_300_3ext, "300_comp/grass_cb_vs_ext", line_data, line_rv, "300", "300+3ext")
# plot_bin(GRASS_rv_300, GRASS_rv_300_3ext, rv_error_GRASS_cb, grass_data_300, grass_data_300_3ext, "300_comp/grass_cb_vs_ext_binned", line_data, line_rv, "300", "300+3ext")
# plot_single_rms(line_data, line_rv, file_regular_300, RV_list_no_cb_300, grass_data_300, GRASS_rv_300, grass_data_no_cb_300, GRASS_no_cb_v_300,
#                 file_regular_300_3ext, RV_list_no_cb_300_3ext, grass_data_300_3ext, GRASS_rv_300_3ext, grass_data_no_cb_300_3ext, GRASS_no_cb_v_300_3ext,
#                 "300_comp/single_rms")

# plot_2_2(GRASS_rv_300_3ext, GRASS_rv_SSD_3ext, rv_error_GRASS_cb, grass_data_300_3ext, grass_data_SSD_3ext, "LD_comp/ext_comp", line_data, line_rv, "300+3ext", "SSD+3ext")

# rms_comp(line_data, line_rv, "RMS_comp")