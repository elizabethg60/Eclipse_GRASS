import h5py
import numpy as np
import pandas as pd
import csv
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from barycorrpy import get_BC_vel

avg_photon_noise = [1.8211758584280178, 1.95054883079038, 2.144446659177623, 2.1658907566938175, 2.004090320949803, 2.596217972425998, 2.1135679819576634, 2.1978020386922146, 2.7521158957968668, 2.1516080460065217, 2.1281340447634953, 2.156445219458067, 3.029994032965276, 2.688447833232973, 3.1333750313001123, 3.341592265251292, 3.9849485492751655, 4.2694904053037765, 3.9731973397147713, 3.6154228169078944, 4.429243144558249, 4.150029318259546]
in_photon_noise = [1.996277941214977, 2.1372306590061374, 2.3488276799304533, 2.372230474620801, 2.194698997502396, 2.840361808394342, 2.3158024928539045, 2.406179964387821, 3.014259096417551, 2.3564576526291905, 2.3304859030659326, 2.3619017742472015, 3.316240627006103, 2.944071902720959, 3.424395677644273, 3.653312371272834, 4.355037342274207, 4.665631965844222, 4.341473323201208, 3.952184601674581, 4.837279462785472, 4.533946285494229]
out_photon_noise = [1.3529152212382396, 1.4513358494944395, 1.5979087739210869, 1.6140965712461186, 1.4943834122812576, 1.943350870898378, 1.5727773284686484, 1.6405873495066787, 2.0511217547762675, 1.6038008223334244, 1.5870162191099135, 1.6070095453998496, 2.2645233808554694, 2.004878299147675, 2.3551159847019085, 2.5080307909702446, 2.995275109346105, 3.2101270772559527, 2.9883342223897453, 2.7148813125244593, 3.3380948809935758, 3.1233459114997477]

def bin_array(arr, bin_size):
    arr = np.array(arr)
    n_bins = len(arr) // bin_size
    return np.array([np.nanmean(arr[i*bin_size:(i+1)*bin_size]) for i in range(n_bins)])

def rms(arr, model):
    return np.sqrt((np.nansum((arr - model)**2))/len(arr))   

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
time_julian = []
UTC_time = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
UTC_time = UTC_time[0:-28]

# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
lines = line_data["name"][()]
line_rv = line_data["rv"][()]

# LD
# projected rv model - regular
file_regular_LD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD/data/projected_SSD_4parameter_gpu.jld2", "r")
RV_list_no_cb_LD = file_regular_LD["RV_list_no_cb"][()]
# GRASS var on
grass_data_LD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD/data/neid_all_lines_rv_on_SSD_4parameter_gpu.jld2", "r")
GRASS_rv_LD  = grass_data_LD["rv"][()]
# GRASS var off
grass_data_no_cb_LD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD/data/neid_all_lines_rv_off_SSD_4parameter_gpu.jld2", "r")
GRASS_no_cb_v_LD  = grass_data_no_cb_LD["rv"][()]

#LD + Ext
# projected rv model - regular
file_regular_LD_ext= h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD_Ext/data/projected_SSD_4parameter_gpu_ext.jld2", "r")
RV_list_no_cb_LD_ext = file_regular_LD_ext["RV_list_no_cb"][()]
# GRASS var on
grass_data_LD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD_Ext/data/neid_all_lines_rv_on_SSD_4parameter_gpu_ext.jld2", "r")
GRASS_rv_LD_ext  = grass_data_LD_ext["rv"][()]
# GRASS var off
grass_data_no_cb_LD_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD_Ext/data/neid_all_lines_rv_off_SSD_4parameter_gpu_ext.jld2", "r")
GRASS_no_cb_v_LD_ext  = grass_data_no_cb_LD_ext["rv"][()]

#Optim CB
# projected rv model - regular
file_regular_LD_CB = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB/data/projected_SSD_4parameter_gpu_ext_CB_optim.jld2", "r")
RV_list_no_cb_LD_CB = file_regular_LD_CB["RV_list_no_cb"][()]
# GRASS var off
grass_data_no_cb_LD_CB = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB/data/neid_all_lines_rv_off_SSD_4parameter_gpu_ext_CB_optim.jld2", "r")
GRASS_no_cb_v_LD_CB  = grass_data_no_cb_LD_CB["rv"][()]

#Optim CB + MF
# projected rv model - regular
file_regular_LD_MF = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB_MF/data/projected_SSD_4parameter_gpu_ext_CB_MF_optim.jld2", "r")
RV_list_no_cb_LD_MF = file_regular_LD_MF["RV_list_no_cb"][()]
# GRASS var off
grass_data_no_cb_LD_MF = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB_MF/data/neid_all_lines_rv_off_SSD_4parameter_gpu_ext_CB_MF_optim.jld2", "r")
GRASS_no_cb_v_LD_MF  = grass_data_no_cb_LD_MF["rv"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_3(projected, grass_cb, grass_no_cb, projected_file, grass_cb_file, grass_no_cb_file, title, label1, label2, label3):
    
    projected_arr = []
    grass_cb_arr = []
    grass_no_cb_arr = []

    projected_arr_binned = []
    grass_cb_arr_binned = []
    grass_no_cb_arr_binned = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_data, line_rv, vb, i)

        projected_arr.append(rms(line_i, projected_i))
        grass_cb_arr.append(rms(line_i, grass_cb_i))
        grass_no_cb_arr.append(rms(line_i, grass_no_cb_i))

        projected_arr_binned.append(rms(bin_array(line_i,4), bin_array(projected_i,4)))
        grass_cb_arr_binned.append(rms(bin_array(line_i,4), bin_array(grass_cb_i,4)))
        grass_no_cb_arr_binned.append(rms(bin_array(line_i,4), bin_array(grass_no_cb_i,4)))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines, projected_arr, label = label1, color = 'b')
    plt.scatter(lines, grass_cb_arr, label = label2, color = 'r')
    plt.scatter(lines, grass_no_cb_arr, label = label3, color = 'g')
    plt.scatter(lines, projected_arr_binned, marker="o", facecolors="none", edgecolors= 'b')
    plt.scatter(lines, grass_cb_arr_binned, marker="o", facecolors="none", edgecolors= 'r')
    plt.scatter(lines, grass_no_cb_arr_binned, marker="o", facecolors="none", edgecolors= 'g')
    plt.scatter(lines, avg_photon_noise, color = 'k', marker='x')
    plt.title("solid - unbinned line rms & circle - binned line rms")
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_3_binned_only(projected, grass_cb, grass_no_cb, projected_file, grass_cb_file, grass_no_cb_file, title, label1, label2, label3):

    projected_arr_binned = []
    grass_cb_arr_binned = []
    grass_no_cb_arr_binned = []
    for i in range(0,len(lines)):

        projected_i = jld2_read(projected_file, projected, vb, i)
        grass_cb_i = jld2_read(grass_cb_file, grass_cb, vb, i)
        grass_no_cb_i = jld2_read(grass_no_cb_file, grass_no_cb, vb, i)
        line_i = jld2_read(line_data, line_rv, vb, i)

        projected_arr_binned.append(rms(bin_array(line_i,4), bin_array(projected_i,4)))
        grass_cb_arr_binned.append(rms(bin_array(line_i,4), bin_array(grass_cb_i,4)))
        grass_no_cb_arr_binned.append(rms(bin_array(line_i,4), bin_array(grass_no_cb_i,4)))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines, projected_arr_binned, marker="o", facecolors="none", edgecolors= 'b', label = label1)
    plt.scatter(lines, grass_cb_arr_binned, marker="o", facecolors="none", edgecolors= 'r', label = label2)
    plt.scatter(lines, grass_no_cb_arr_binned, marker="o", facecolors="none", edgecolors= 'g', label = label3)
    plt.scatter(lines, avg_photon_noise, color = 'k', marker='x', label = "Photon Noise")
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(title, bbox_inches='tight')
    plt.clf() 

def plot_2(projected1, grass_no_cb1, projected_file1, grass_no_cb_file1, projected2, grass_no_cb2, projected_file2, grass_no_cb_file2, projected3, grass_no_cb3, projected_file3, grass_no_cb_file3, title, savefig, label1, label2):
    
    projected_arr1 = []
    grass_no_cb_arr1 = []

    projected_arr2 = []
    grass_no_cb_arr2 = []

    projected_arr3 = []
    grass_no_cb_arr3 = []
    for i in range(0,len(lines)):

        projected_i1 = jld2_read(projected_file1, projected1, vb, i)
        grass_no_cb_i1 = jld2_read(grass_no_cb_file1, grass_no_cb1, vb, i)
        projected_i2 = jld2_read(projected_file2, projected2, vb, i)
        grass_no_cb_i2 = jld2_read(grass_no_cb_file2, grass_no_cb2, vb, i)
        projected_i3 = jld2_read(projected_file3, projected3, vb, i)
        grass_no_cb_i3 = jld2_read(grass_no_cb_file3, grass_no_cb3, vb, i)
        line_i = jld2_read(line_data, line_rv, vb, i)

        projected_arr1.append(rms(line_i, projected_i1))
        grass_no_cb_arr1.append(rms(line_i, grass_no_cb_i1))

        projected_arr2.append(rms(line_i, projected_i2))
        grass_no_cb_arr2.append(rms(line_i, grass_no_cb_i2))

        projected_arr3.append(rms(line_i, projected_i3))
        grass_no_cb_arr3.append(rms(line_i, grass_no_cb_i3))

    plt.figure(figsize=(12, 6))
    plt.scatter(lines, projected_arr1, label = label1, color = 'b')
    plt.scatter(lines, grass_no_cb_arr1, label = label2, color = 'g')
    plt.scatter(lines, projected_arr2, marker="o", facecolors="none", edgecolors= 'b')
    plt.scatter(lines, grass_no_cb_arr2, marker="o", facecolors="none", edgecolors= 'g')
    plt.scatter(lines, projected_arr3, marker="s", color = 'b')
    plt.scatter(lines, grass_no_cb_arr3, marker="s", color = 'g')
    plt.scatter(lines, avg_photon_noise, color = 'k', marker='x')
    plt.title(title)
    plt.xlabel("Line Wavelength (Å)", fontsize=12)
    plt.ylabel("RV RMS (m/s)", fontsize=12)
    plt.xticks(rotation=60)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(savefig, bbox_inches='tight')
    plt.clf() 

# plot_3(RV_list_no_cb_LD, GRASS_rv_LD, GRASS_no_cb_v_LD, file_regular_LD, grass_data_LD, grass_data_no_cb_LD, "LD_comp", "model one", "model two", "model three")
# plot_3(RV_list_no_cb_LD_ext, GRASS_rv_LD_ext, GRASS_no_cb_v_LD_ext, file_regular_LD_ext, grass_data_LD_ext, grass_data_no_cb_LD_ext, "Ext_comp", "model one", "model two", "model three")
# plot_3_binned_only(RV_list_no_cb_LD_CB, GRASS_rv_LD_ext, GRASS_no_cb_v_LD_CB, file_regular_LD_CB, grass_data_LD_ext, grass_data_no_cb_LD_CB, "CB_comps_binned_only", "Model One", "Model Two", "Model Three")
# plot_3(RV_list_no_cb_LD_CB, GRASS_rv_LD_ext, GRASS_no_cb_v_LD_CB, file_regular_LD_CB, grass_data_LD_ext, grass_data_no_cb_LD_CB, "CB_comps", "model one (with CB)", "model two", "model three (with CB)")
plot_3_binned_only(RV_list_no_cb_LD_MF, GRASS_rv_LD_ext, GRASS_no_cb_v_LD_MF, file_regular_LD_MF, grass_data_LD_ext, grass_data_no_cb_LD_MF, "CB_MF_comps_binned_only", "Model One", "Model Two", "Model Three")
plot_3(RV_list_no_cb_LD_MF, GRASS_rv_LD_ext, GRASS_no_cb_v_LD_MF, file_regular_LD_MF, grass_data_LD_ext, grass_data_no_cb_LD_MF, "CB_MF_comps", "model one (with CB & MF)", "model two", "model three (with CB & MF)")
plot_2(RV_list_no_cb_LD_ext, GRASS_no_cb_v_LD_ext, file_regular_LD_ext, grass_data_no_cb_LD_ext, RV_list_no_cb_LD_CB, GRASS_no_cb_v_LD_CB, file_regular_LD_CB, grass_data_no_cb_LD_CB, 
        RV_list_no_cb_LD_MF, GRASS_no_cb_v_LD_MF, file_regular_LD_MF, grass_data_no_cb_LD_MF,
        "solid - no CB rms & circle - optimized CB rms & square - optimized CB and MF rms", "Optim_comps", "model one", "model three")