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

def rms(arr, model):
    return np.sqrt((np.nansum((arr - model)**2))/len(arr))    

def weighted_rms(arr, model, weight):
    # weight = []
    # for i in range(0, len(arr)):
    #     weight.append(np.sqrt((np.nansum((arr[i] - model[i])**2))))
    return np.sqrt((np.nansum((arr - model)**2 * weight))/(np.nansum(weight)*len(arr)))

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
time_julian = []
UTC_time = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    #UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S:%f"))
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)[0:-25]

df = pd.read_csv("/storage/home/efg5335/work/GRASS/data/optimized_depth.csv")
df = df.drop(10)
df = df.reset_index(drop=True)


#NL94
#projected rv model - regular
file_regular_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_october_N_50_NL94.jld2", "r")
RV_list_no_cb_NL94 = file_regular_NL94["RV_list_no_cb"][()]
#GRASS CB
grass_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_all_lines_rv_regular_NL94.jld2", "r")
GRASS_rv_NL94  = grass_data_NL94["rv"][()]
#GRASS no CB
grass_data_no_cb_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/NL94/data/neid_all_lines_rv_off_NL94.jld2", "r")
GRASS_no_cb_v_NL94  = grass_data_no_cb_NL94["rv"][()]
#line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
lines = line_data["name"][()]
line_rv = line_data["rv"][()]
line_rv_err = line_data["rv_error"][()]

#K300
file_regular_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300/data/neid_october_N_50_K300.jld2", "r")
RV_list_no_cb_300 = file_regular_300["RV_list_no_cb"][()]
#GRASS CB
grass_data_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300/data/neid_all_lines_rv_regular_K300.jld2", "r")
GRASS_rv_300  = grass_data_300["rv"][()]
#GRASS no CB
grass_data_no_cb_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300/data/neid_all_lines_rv_off_K300.jld2", "r")
GRASS_no_cb_v_300  = grass_data_no_cb_300["rv"][()]

#K300 + 3ext
file_regular_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300_3ext/data/neid_october_N_50_K300_3ext.jld2", "r")
RV_list_no_cb_300_3ext = file_regular_300_3ext["RV_list_no_cb"][()]
#GRASS CB
grass_data_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300_3ext/data/neid_all_lines_rv_regular_K300_3ext.jld2", "r")
GRASS_rv_300_3ext  = grass_data_300_3ext["rv"][()]
#GRASS no CB
grass_data_no_cb_300_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/K300/K300_3ext/data/neid_all_lines_rv_off_K300_3ext.jld2", "r")
GRASS_no_cb_v_300_3ext  = grass_data_no_cb_300_3ext["rv"][()]

#KSSD
#projected rv model - regular
file_regular_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD/data/neid_october_N_50_KSSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_regular_KSSD.jld2", "r")
GRASS_rv_KSSD  = grass_data_KSSD["rv"][()]
rv_error_GRASS_cb  = grass_data_KSSD["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_off_KSSD.jld2", "r")
GRASS_no_cb_v_KSSD  = grass_data_no_cb_KSSD["rv"][()]

#KSSD with 3 extinction 
#projected rv model - regular
file_regular_SSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD_3ext/data/neid_october_N_50_KSSD_3ext.jld2", "r")
RV_list_no_cb_SSD_3ext = file_regular_SSD_3ext["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD_3ext/data/neid_all_lines_rv_regular_KSSD_3ext.jld2", "r")
GRASS_rv_KSSD_3ext  = grass_data_KSSD_3ext["rv"][()]
rv_error_GRASS_cb_3ext  = grass_data_KSSD_3ext["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD_3ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD_3ext/data/neid_all_lines_rv_off_KSSD_3ext.jld2", "r")
GRASS_no_cb_v_KSSD_3ext  = grass_data_no_cb_KSSD_3ext["rv"][()]

#HD
#projected rv model - regular
file_regular_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KHD/data/neid_october_N_50_KHD.jld2", "r")
RV_list_no_cb_HD = file_regular_HD["RV_list_no_cb"][()]
#GRASS CB
grass_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KHD/data/neid_all_lines_rv_regular_KHD.jld2", "r")
GRASS_rv_HD  = grass_data_HD["rv"][()]
#GRASS no CB
grass_data_no_cb_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KHD/data/neid_all_lines_rv_off_KHD.jld2", "r")
GRASS_no_cb_v_HD  = grass_data_no_cb_HD["rv"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    array -= array[-1]
    return array

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
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-25]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))

        projected_arr_weight.append(weighted_rms(line_i, projected_i, line_file[line_rv_err[i]][()][0:-25]))
        grass_cb_arr_weight.append(weighted_rms(line_i, grass_cb_i, line_file[line_rv_err[i]][()][0:-25]))
        grass_no_cb_arr_weight.append(weighted_rms(line_i, grass_no_cb_i, line_file[line_rv_err[i]][()][0:-25]))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))
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
    if title == "KSSD" or title == "grass_cb_vs_ld":
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
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        data_300_i = jld2_read(file_300, data_300, vb, i)
        
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-25]))

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))
        arr_300.append(rms(line_i, data_300_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))
        arr_300_pipe.append(rms(rv_obs, data_300_i))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))

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

    if title == "KSSD" or title == "grass_cb_vs_ld":
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
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-25]))

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
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)
        line_err.append(np.mean(line_file[line_rv_err[i]][()][0:-25]))

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
    KSSD_projected = []
    KSSD_grass_cb = []
    KSSD_grass_no_cb = []
    KSSD_3ext_projected = []
    KSSD_3ext_grass_cb = []
    KSSD_3ext_grass_no_cb = []
    KSSD_projected_pipeline = []
    KSSD_grass_cb_pipeline = []
    KSSD_grass_no_cb_pipeline = []
    KSSD_3ext_projected_pipeline = []
    KSSD_3ext_grass_cb_pipeline = []
    KSSD_3ext_grass_no_cb_pipeline = []

    for i in range(0,len(lines)):
        if i == 10:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        projected_i_3ext = jld2_read(projected_file_3ext, projected_3ext, vb, i)
        grass_cb_i_3ext = jld2_read(grass_cb_file_3ext, grass_cb_3ext, vb, i)
        grass_no_cb_i_3ext = jld2_read(grass_no_cb_file_3ext, grass_no_cb_3ext, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)

        for j in range(0, len(line_i)):
            KSSD_projected.append(line_i[j] - projected_i[j])
            KSSD_grass_cb.append(line_i[j] - grass_cb_i[j])
            KSSD_grass_no_cb.append(line_i[j] - grass_no_cb_i[j])
            KSSD_3ext_projected.append(line_i[j] - projected_i_3ext[j])
            KSSD_3ext_grass_cb.append(line_i[j] - grass_cb_i_3ext[j])
            KSSD_3ext_grass_no_cb.append(line_i[j] - grass_no_cb_i_3ext[j])

            KSSD_projected_pipeline.append(rv_obs[j] - projected_i[j])
            KSSD_grass_cb_pipeline.append(rv_obs[j] - grass_cb_i[j])
            KSSD_grass_no_cb_pipeline.append(rv_obs[j] - grass_no_cb_i[j])
            KSSD_3ext_projected_pipeline.append(rv_obs[j] - projected_i_3ext[j])
            KSSD_3ext_grass_cb_pipeline.append(rv_obs[j] - grass_cb_i_3ext[j])
            KSSD_3ext_grass_no_cb_pipeline.append(rv_obs[j] - grass_no_cb_i_3ext[j])

    KSSD_projected = np.array(KSSD_projected)
    KSSD_grass_cb = np.array(KSSD_grass_cb)
    KSSD_grass_no_cb = np.array(KSSD_grass_no_cb)
    KSSD_3ext_projected = np.array(KSSD_3ext_projected)
    KSSD_3ext_grass_cb = np.array(KSSD_3ext_grass_cb)
    KSSD_3ext_grass_no_cb = np.array(KSSD_3ext_grass_no_cb)

    KSSD_projected_pipeline = np.array(KSSD_projected_pipeline)
    KSSD_grass_cb_pipeline = np.array(KSSD_grass_cb_pipeline)
    KSSD_grass_no_cb_pipeline = np.array(KSSD_grass_no_cb_pipeline)
    KSSD_3ext_projected_pipeline = np.array(KSSD_3ext_projected_pipeline)
    KSSD_3ext_grass_cb_pipeline = np.array(KSSD_3ext_grass_cb_pipeline)
    KSSD_3ext_grass_no_cb_pipeline = np.array(KSSD_3ext_grass_no_cb_pipeline)

    x = [0,1,2]   
    y_kssd = [np.sqrt((np.nansum((KSSD_projected)**2))/len(KSSD_projected)), np.sqrt((np.nansum((KSSD_grass_cb)**2))/len(KSSD_grass_cb)), np.sqrt((np.nansum((KSSD_grass_no_cb)**2))/len(KSSD_grass_no_cb))]
    y_kssd_3ext = [np.sqrt((np.nansum((KSSD_3ext_projected)**2))/len(KSSD_3ext_projected)), np.sqrt((np.nansum((KSSD_3ext_grass_cb)**2))/len(KSSD_3ext_grass_cb)), np.sqrt((np.nansum((KSSD_3ext_grass_no_cb)**2))/len(KSSD_3ext_grass_no_cb))]
    y_kssd_pipeline = [np.sqrt((np.nansum((KSSD_projected_pipeline)**2))/len(KSSD_projected_pipeline)), np.sqrt((np.nansum((KSSD_grass_cb_pipeline)**2))/len(KSSD_grass_cb_pipeline)), np.sqrt((np.nansum((KSSD_grass_no_cb_pipeline)**2))/len(KSSD_grass_no_cb_pipeline))]
    y_kssd_3ext_pipeline = [np.sqrt((np.nansum((KSSD_3ext_projected_pipeline)**2))/len(KSSD_3ext_projected_pipeline)), np.sqrt((np.nansum((KSSD_3ext_grass_cb_pipeline)**2))/len(KSSD_3ext_grass_cb_pipeline)), np.sqrt((np.nansum((KSSD_3ext_grass_no_cb_pipeline)**2))/len(KSSD_3ext_grass_no_cb_pipeline))]
    labels = ['Projected', 'GRASS Var', 'GRASS No Var']
    plt.figure(figsize=(12, 6))
    plt.scatter(x, y_kssd, color = 'b', label = "KSSD")
    plt.scatter(x, y_kssd_3ext, color = 'g', label = "KSSD+3ext")

    plt.scatter(x, y_kssd_pipeline,  color = 'b', facecolors='none')
    plt.scatter(x, y_kssd_3ext_pipeline, color = 'g', facecolors='none')
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


plot_3(RV_list_no_cb_NL94, GRASS_rv_NL94, GRASS_no_cb_v_NL94, 0, file_regular_NL94, grass_data_NL94, grass_data_no_cb_NL94, "LD_comp/NL94", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_3(RV_list_no_cb_SSD, GRASS_rv_KSSD, GRASS_no_cb_v_KSSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_KSSD, grass_data_no_cb_KSSD, "LD_comp/KSSD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_3(RV_list_no_cb_HD, GRASS_rv_HD, GRASS_no_cb_v_HD, 0, file_regular_HD, grass_data_HD, grass_data_no_cb_HD, "LD_comp/HD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_3(RV_list_no_cb_300, GRASS_rv_300, GRASS_no_cb_v_300, 0, file_regular_300, grass_data_300, grass_data_no_cb_300, "LD_comp/300", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_4_4(GRASS_rv_NL94, GRASS_rv_KSSD, GRASS_rv_HD, GRASS_rv_300, rv_error_GRASS_cb, grass_data_NL94, grass_data_KSSD, grass_data_HD, grass_data_300, "LD_comp/grass_cb_vs_ld", line_data, line_rv, "NL94", "SSD", "HD", "300")

plot_3(RV_list_no_cb_SSD, GRASS_rv_KSSD, GRASS_no_cb_v_KSSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_KSSD, grass_data_no_cb_KSSD, "KSSD_comp/KSSD", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_3(RV_list_no_cb_SSD_3ext, GRASS_rv_KSSD_3ext, GRASS_no_cb_v_KSSD_3ext, rv_error_GRASS_cb_3ext, file_regular_SSD_3ext, grass_data_KSSD_3ext, grass_data_no_cb_KSSD_3ext, "KSSD_comp/KSSD_3ext", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_2_2(GRASS_rv_KSSD, GRASS_rv_KSSD_3ext, rv_error_GRASS_cb, grass_data_KSSD, grass_data_KSSD_3ext, "KSSD_comp/grass_cb_vs_ext", line_data, line_rv, "SSD", "SSD+3ext")
plot_bin(GRASS_rv_KSSD, GRASS_rv_KSSD_3ext, rv_error_GRASS_cb, grass_data_KSSD, grass_data_KSSD_3ext, "KSSD_comp/grass_cb_vs_ext_binned", line_data, line_rv, "SSD", "SSD+3ext")
plot_single_rms(line_data, line_rv, file_regular_SSD, RV_list_no_cb_SSD, grass_data_KSSD, GRASS_rv_KSSD, grass_data_no_cb_KSSD, GRASS_no_cb_v_KSSD,
                file_regular_SSD_3ext, RV_list_no_cb_SSD_3ext, grass_data_KSSD_3ext, GRASS_rv_KSSD_3ext, grass_data_no_cb_KSSD_3ext, GRASS_no_cb_v_KSSD_3ext,
                "KSSD_comp/single_rms")

plot_3(RV_list_no_cb_300, GRASS_rv_300, GRASS_no_cb_v_300, rv_error_GRASS_cb, file_regular_300, grass_data_300, grass_data_no_cb_300, "K300_comp/K300", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_3(RV_list_no_cb_300_3ext, GRASS_rv_300_3ext, GRASS_no_cb_v_300_3ext, rv_error_GRASS_cb_3ext, file_regular_300_3ext, grass_data_300_3ext, grass_data_no_cb_300_3ext, "K300_comp/K300_3ext", line_data, line_rv, 'Projected RV - no cb', 'CCF RV - GRASS var', 'CCF RV - no var')
plot_2_2(GRASS_rv_300, GRASS_rv_300_3ext, rv_error_GRASS_cb, grass_data_300, grass_data_300_3ext, "K300_comp/grass_cb_vs_ext", line_data, line_rv, "300", "300+3ext")
plot_bin(GRASS_rv_300, GRASS_rv_300_3ext, rv_error_GRASS_cb, grass_data_300, grass_data_300_3ext, "K300_comp/grass_cb_vs_ext_binned", line_data, line_rv, "300", "300+3ext")
plot_single_rms(line_data, line_rv, file_regular_300, RV_list_no_cb_300, grass_data_300, GRASS_rv_300, grass_data_no_cb_300, GRASS_no_cb_v_300,
                file_regular_300_3ext, RV_list_no_cb_300_3ext, grass_data_300_3ext, GRASS_rv_300_3ext, grass_data_no_cb_300_3ext, GRASS_no_cb_v_300_3ext,
                "K300_comp/single_rms")

plot_2_2(GRASS_rv_300, GRASS_rv_KSSD, rv_error_GRASS_cb, grass_data_300, grass_data_KSSD, "LD_comp/no_ext_comp", line_data, line_rv, "300", "SSD")
plot_2_2(GRASS_rv_300_3ext, GRASS_rv_KSSD_3ext, rv_error_GRASS_cb, grass_data_300_3ext, grass_data_KSSD_3ext, "LD_comp/ext_comp", line_data, line_rv, "300+3ext", "SSD+3ext")
