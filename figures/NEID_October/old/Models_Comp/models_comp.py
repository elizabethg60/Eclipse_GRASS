import h5py
import numpy as np
import pandas as pd
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

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)[0:-25]

#NL94
#projected rv model - regular
file_regular_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_october_N_50_NL94.jld2", "r")
RV_list_no_cb_NL94 = file_regular_NL94["RV_list_no_cb"][()]
#GRASS CB
grass_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_regular_NL94.jld2", "r")
GRASS_rv_NL94  = grass_data_NL94["rv"][()]
#GRASS no CB
grass_data_no_cb_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_off_NL94.jld2", "r")
GRASS_no_cb_v_NL94  = grass_data_no_cb_NL94["rv"][()]
#line by line data
line_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_RVlinebyline_NL94.jld2", "r")
lines = line_data_NL94["name"][()]
line_rv_NL94  = line_data_NL94["rv"][()]

#KSSD
#projected rv model - regular
file_regular_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_october_N_50_KSSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_regular_KSSD.jld2", "r")
GRASS_rv_KSSD  = grass_data_KSSD["rv"][()]
rv_error_GRASS_cb  = grass_data_KSSD["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_off_KSSD.jld2", "r")
GRASS_no_cb_v_KSSD  = grass_data_no_cb_KSSD["rv"][()]
#line by line data
line_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_RVlinebyline_KSSD.jld2", "r")
line_rv_KSSD  = line_data_KSSD["rv"][()]

#KSSD with extinction 
#projected rv model - regular
file_regular_SSD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_october_N_50_KSSD_ext.jld2", "r")
RV_list_no_cb_SSD_ext = file_regular_SSD_ext["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_all_lines_rv_regular_KSSD_ext.jld2", "r")
GRASS_rv_KSSD_ext  = grass_data_KSSD_ext["rv"][()]
rv_error_GRASS_cb_ext  = grass_data_KSSD_ext["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_all_lines_rv_off_KSSD_ext.jld2", "r")
GRASS_no_cb_v_KSSD_ext  = grass_data_no_cb_KSSD_ext["rv"][()]
#line by line data
line_data_KSSD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_RVlinebyline_KSSD_ext.jld2", "r")
line_rv_KSSD_ext  = line_data_KSSD_ext["rv"][()]

#KSSD with 2 extinction 
#projected rv model - regular
file_regular_SSD_2ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/data/neid_october_N_50_KSSD_2_ext.jld2", "r")
RV_list_no_cb_SSD_2ext = file_regular_SSD_2ext["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD_2ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/data/neid_all_lines_rv_regular_KSSD_2_ext.jld2", "r")
GRASS_rv_KSSD_2ext  = grass_data_KSSD_2ext["rv"][()]
rv_error_GRASS_cb_2ext  = grass_data_KSSD_2ext["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD_2ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/data/neid_all_lines_rv_off_KSSD_2_ext.jld2", "r")
GRASS_no_cb_v_KSSD_2ext  = grass_data_no_cb_KSSD_2ext["rv"][()]
#line by line data
line_data_KSSD_2ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/data/neid_RVlinebyline_KSSD_2_ext.jld2", "r")
line_rv_KSSD_2ext  = line_data_KSSD_2ext["rv"][()]

#KSSD with 2 extinction + geometry changes
#projected rv model - regular
file_regular_SSD_2ext_new = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/Reiners_Changes/data/neid_october_N_50_KSSD_2_ext_changes2.jld2", "r")
RV_list_no_cb_SSD_2ext_new = file_regular_SSD_2ext_new["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD_2ext_new = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/Reiners_Changes/data/neid_all_lines_rv_regular_KSSD_2_ext_changes2.jld2", "r")
GRASS_rv_KSSD_2ext_new  = grass_data_KSSD_2ext_new["rv"][()]
rv_error_GRASS_cb_2ext_new  = grass_data_KSSD_2ext_new["rv_error"][()]
#GRASS no CB
grass_data_no_cb_KSSD_2ext_new = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/Reiners_Changes/data/neid_all_lines_rv_off_KSSD_2_ext_changes2.jld2", "r")
GRASS_no_cb_v_KSSD_2ext_new  = grass_data_no_cb_KSSD_2ext_new["rv"][()]
#line by line data
line_data_KSSD_2ext_new = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_2ext/Reiners_Changes/data/neid_RVlinebyline_KSSD_2_ext_changes2.jld2", "r")
line_rv_KSSD_2ext_new  = line_data_KSSD_2ext_new["rv"][()]

#HD
#projected rv model - regular
file_regular_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_october_N_50_HD.jld2", "r")
RV_list_no_cb_HD = file_regular_HD["RV_list_no_cb"][()]
#GRASS CB
grass_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_all_lines_rv_regular_HD.jld2", "r")
GRASS_rv_HD  = grass_data_HD["rv"][()]
#GRASS no CB
grass_data_no_cb_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_all_lines_rv_off_HD.jld2", "r")
GRASS_no_cb_v_HD  = grass_data_no_cb_HD["rv"][()]
#line by line data
line_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_RVlinebyline_HD.jld2", "r")
line_rv_HD  = line_data_HD["rv"][()]

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
    lines_arr = []
    grass_cb_err_arr = []
    for i in range(0,len(lines)):
        if i == 5:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_pipe, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none')
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

def plot_3_3(projected, grass_cb, grass_no_cb, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, title, 
                 line_file1, line_rv1, line_file2, line_rv2, line_file3, line_rv3,
                 label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    projected_arr_pipe = []
    grass_cb_arr_pipe = []
    grass_no_cb_arr_pipe = []
    lines_arr = []
    grass_cb_err_arr = []
    for i in range(0,len(lines)):
        if i == 5:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i1 = jld2_read(line_file1, line_rv1, vb, i)
        line_i2 = jld2_read(line_file2, line_rv2, vb, i)
        line_i3 = jld2_read(line_file3, line_rv3, vb, i)

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i1, projected_i))
        grass_cb_arr.append(rms(line_i2, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i3, grass_no_cb_i))

        projected_arr_pipe.append(rms(rv_obs, projected_i))
        grass_cb_arr_pipe.append(rms(rv_obs, grass_cb_i))
        grass_no_cb_arr_pipe.append(rms(rv_obs, grass_no_cb_i))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_pipe, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_pipe, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_pipe, color = 'g', facecolors='none')
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

def plot_bin(projected, grass_cb, grass_no_cb, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, title, line_file1, line_rv1,
                 line_file2, line_rv2, line_file3, line_rv3,
                 label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    projected_arr_bin = []
    grass_cb_arr_bin = []
    grass_no_cb_arr_bin = []
    lines_arr = []
    grass_cb_err_arr = []
    for i in range(0,len(lines)):
        if i == 5:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i1 = jld2_read(line_file1, line_rv1, vb, i)
        line_i2 = jld2_read(line_file2, line_rv2, vb, i)
        line_i3 = jld2_read(line_file3, line_rv3, vb, i)

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i1, projected_i))
        grass_cb_arr.append(rms(line_i2, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i3, grass_no_cb_i))

        projected_arr_bin.append(rms(bin_array(line_i1, 5), bin_array(projected_i,5)))
        grass_cb_arr_bin.append(rms(bin_array(line_i2, 5), bin_array(grass_cb_i,5)))
        grass_no_cb_arr_bin.append(rms(bin_array(line_i3, 5), bin_array(grass_no_cb_i,5)))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_bin, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_bin, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_bin, color = 'g', facecolors='none')
    if title == "KSSD" or title == "grass_cb_vs_ld":
        plt.scatter(lines_arr, grass_cb_err_arr, color = 'k', marker = 'x')
        plt.title("solid - line rms & circle - bin rms & x - mean error grass cb ccf (SSD)")
    else:
        plt.title("solid - line rms & circle - bin rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf()     

def plot_6(projected, grass_cb, grass_no_cb, grass_cb_err, 
                 projected_file, grass_cb_file, grass_no_cb_file, line_file, line_rv,
                 projected_new, grass_cb_new, grass_no_cb_new, grass_cb_err_new, 
                 projected_file_new, grass_cb_file_new, grass_no_cb_file_new, title, line_file_new, line_rv_new,
                 label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []
    projected_arr_new = []
    grass_cb_arr_new = []
    grass_no_cb_arr_new = []
    lines_arr = []
    grass_cb_err_arr = []
    grass_cb_err_arr_new = []

    for i in range(0,len(lines)):
        if i == 5:
            continue

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_file, line_rv, vb, i)

        projected_i_new = jld2_read(projected_file_new, projected_new, vb, i)
        grass_cb_i_new = jld2_read(grass_cb_file_new, grass_cb_new, vb, i)
        grass_no_cb_i_new = jld2_read(grass_no_cb_file_new, grass_no_cb_new, vb, i)
        line_i_new = jld2_read(line_file_new, line_rv_new, vb, i)

        lines_arr.append(lines[i])
        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))
        projected_arr_new.append(rms(line_i_new, projected_i_new))
        grass_cb_arr_new.append(rms(line_i_new, grass_cb_i_new))
        grass_no_cb_arr_new.append(rms(line_i_new, grass_no_cb_i_new))

        if title == "KSSD" or title == "grass_cb_vs_ld":
            grass_cb_err_arr.append(np.mean(grass_cb_file[grass_cb_err[i]][()][0:-25]))
            grass_cb_err_arr_new.append(np.mean(grass_cb_file_new[grass_cb_err_new[i]][()][0:-25]))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines_arr, projected_arr, label = label1, color = 'b')
    plt.scatter(lines_arr, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines_arr, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines_arr, projected_arr_new, color = 'b', facecolors='none')
    plt.scatter(lines_arr, grass_cb_arr_new, color = 'r', facecolors='none')
    plt.scatter(lines_arr, grass_no_cb_arr_new, color = 'g', facecolors='none')
    if title == "KSSD" or title == "grass_cb_vs_ld":
        plt.scatter(lines_arr, grass_cb_err_arr, color = 'k', marker = 'x')
        plt.scatter(lines_arr, grass_cb_err_arr_new, color = 'grey', marker = 'x')

        plt.title("solid - line rms & circle - changes rms & x - mean error grass cb ccf (SSD)")
    else:
        plt.title("solid - line rms & circle - changes rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf()

plot_3(RV_list_no_cb_NL94, GRASS_rv_NL94, GRASS_no_cb_v_NL94, 0, file_regular_NL94, grass_data_NL94, grass_data_no_cb_NL94, "LD_comp/NL94", line_data_NL94, line_rv_NL94, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3(RV_list_no_cb_SSD, GRASS_rv_KSSD, GRASS_no_cb_v_KSSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_KSSD, grass_data_no_cb_KSSD, "LD_comp/KSSD", line_data_KSSD, line_rv_KSSD, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3(RV_list_no_cb_HD, GRASS_rv_HD, GRASS_no_cb_v_HD, 0, file_regular_HD, grass_data_HD, grass_data_no_cb_HD, "LD_comp/HD", line_data_HD, line_rv_HD, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3_3(GRASS_rv_NL94, GRASS_rv_KSSD, GRASS_rv_HD, rv_error_GRASS_cb, grass_data_NL94, grass_data_KSSD, grass_data_HD, "LD_comp/grass_cb_vs_ld", line_data_NL94, line_rv_NL94, line_data_KSSD, line_rv_KSSD, line_data_HD, line_rv_HD, "NL94", "SSD", "HD")

plot_3(RV_list_no_cb_SSD, GRASS_rv_KSSD, GRASS_no_cb_v_KSSD, rv_error_GRASS_cb, file_regular_SSD, grass_data_KSSD, grass_data_no_cb_KSSD, "KSSD_comp/KSSD", line_data_KSSD, line_rv_KSSD, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3(RV_list_no_cb_SSD_ext, GRASS_rv_KSSD_ext, GRASS_no_cb_v_KSSD_ext, rv_error_GRASS_cb_ext, file_regular_SSD_ext, grass_data_KSSD_ext, grass_data_no_cb_KSSD_ext, "KSSD_comp/KSSD_ext", line_data_KSSD_ext, line_rv_KSSD_ext, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3(RV_list_no_cb_SSD_2ext, GRASS_rv_KSSD_2ext, GRASS_no_cb_v_KSSD_2ext, rv_error_GRASS_cb_2ext, file_regular_SSD_2ext, grass_data_KSSD_2ext, grass_data_no_cb_KSSD_2ext, "KSSD_comp/KSSD_2ext", line_data_KSSD_2ext, line_rv_KSSD_2ext, 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')
plot_3_3(GRASS_rv_KSSD, GRASS_rv_KSSD_ext, GRASS_rv_KSSD_2ext, rv_error_GRASS_cb, grass_data_KSSD, grass_data_KSSD_ext, grass_data_KSSD_2ext, "KSSD_comp/grass_cb_vs_ext", line_data_KSSD, line_rv_KSSD, line_data_KSSD_ext, line_rv_KSSD_ext,line_data_KSSD_2ext, line_rv_KSSD_2ext, "SSD", "SSD+ext", "SSD+2ext")

plot_bin(GRASS_rv_KSSD, GRASS_rv_KSSD_ext, GRASS_rv_KSSD_2ext, rv_error_GRASS_cb, grass_data_KSSD, grass_data_KSSD_ext, grass_data_KSSD_2ext, "KSSD_comp/grass_cb_vs_ext_binned", line_data_KSSD, line_rv_KSSD, line_data_KSSD_ext, line_rv_KSSD_ext,line_data_KSSD_2ext, line_rv_KSSD_2ext, "SSD", "SSD+ext", "SSD+2ext")

plot_6(RV_list_no_cb_SSD_2ext, GRASS_rv_KSSD_2ext, GRASS_no_cb_v_KSSD_2ext, rv_error_GRASS_cb_2ext, 
                 file_regular_SSD_2ext, grass_data_KSSD_2ext, grass_data_no_cb_KSSD_2ext, line_data_KSSD_2ext, line_rv_KSSD_2ext,
                 RV_list_no_cb_SSD_2ext_new, GRASS_rv_KSSD_2ext_new, GRASS_no_cb_v_KSSD_2ext_new, rv_error_GRASS_cb_2ext_new, 
                 file_regular_SSD_2ext_new, grass_data_KSSD_2ext_new, grass_data_no_cb_KSSD_2ext_new,  
                 "KSSD_comp/grass_cb_vs_2ext_changes", line_data_KSSD_2ext_new, line_rv_KSSD_2ext_new,
                 'Projected RV - no cb', 'CCF RV - GRASS cb', 'CCF RV - no cb')