import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

# read in data
# GRASS
grass_data = h5py.File("data/neid_all_lines_rv_regular_SSD_gpu_first.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
rv_error_GRASS_cb  = grass_data["rv_error"][()]
grass_data_no_cb = h5py.File("data/neid_all_lines_rv_off_SSD_gpu_first.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
rv_error_GRASS_no_cb = grass_data_no_cb["rv_error"][()]

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
rv_error_line = line_data["rv_error"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs[0:-28])
rv_obs -= rv_obs[-1]

UTC_time = UTC_time[0:-28]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-28]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_line(UTC_time, rv_obs, line_rv_array, GRASS, GRASS_label, path, save, rv_error_GRASS_no_cb_array, rv_error_line_array):
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID line RVs") 
    axs[0].plot(UTC_time, GRASS, color = 'b', linewidth = 2, label = GRASS_label)
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS)**2))/len(line_rv_array - GRASS)),2)
    axs[0].text(UTC_time[-40], 700, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, line_rv_array - GRASS, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("{}/rm_and_residuals_{}.png".format(path, save))
    plt.clf()

for i in range(0,11):
    rv_error_GRASS_cb_array = grass_data[rv_error_GRASS_cb[i]][()][0:-28]
    rv_error_GRASS_no_cb_array = grass_data_no_cb[rv_error_GRASS_no_cb[i]][()][0:-28]
    rv_error_line_array = line_data[rv_error_line[i]][()][0:-28]
    
    GRASS_rv_array = jld2_read(grass_data, GRASS_rv, vb, i)
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    GRASS_no_cb_array = jld2_read(grass_data_no_cb, GRASS_no_cb, vb, i)

    # rm curve 
    plot_line(UTC_time, rv_obs, line_rv_array, GRASS_no_cb_array, "Line RVs - No Var", "No_Granulation/LinebyLine", lines_no_cb[i], rv_error_GRASS_no_cb_array, rv_error_line_array)

    # rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No Var")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'b', marker = "x", s = 3) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("No_Granulation/rm_and_residuals_{}.png".format(lines_no_cb[i]))
    plt.clf()

    # rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID line RVs")
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_rv_array)**2))/len(line_rv_array - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-40], -300, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals 
    axs[1].scatter(UTC_time, line_rv_array - GRASS_rv_array, color = 'b', marker = "x", s = 3) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("Granulation/LinebyLine/rm_and_residuals_cb_{}.png".format(lines[i]))
    plt.clf()

    # rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)    
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("Granulation/rm_and_residuals_cb_{}.png".format(lines[i]))
    plt.clf()

grass_data = h5py.File("data/neid_all_lines_rv_regular_SSD_gpu.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
rv_error_GRASS_cb  = grass_data["rv_error"][()]
grass_data_no_cb = h5py.File("data/neid_all_lines_rv_off_SSD_gpu.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
rv_error_GRASS_no_cb = grass_data_no_cb["rv_error"][()]

for i in range(11,22):
    rv_error_GRASS_cb_array = grass_data[rv_error_GRASS_cb[i]][()][0:-28]
    rv_error_GRASS_no_cb_array = grass_data_no_cb[rv_error_GRASS_no_cb[i]][()][0:-28]
    rv_error_line_array = line_data[rv_error_line[i]][()][0:-28]
    
    GRASS_rv_array = jld2_read(grass_data, GRASS_rv, vb, i)
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    GRASS_no_cb_array = jld2_read(grass_data_no_cb, GRASS_no_cb, vb, i)

    # rm curve 
    plot_line(UTC_time, rv_obs, line_rv_array, GRASS_no_cb_array, "Line RVs - No Var", "No_Granulation/LinebyLine", lines_no_cb[i], rv_error_GRASS_no_cb_array, rv_error_line_array)

    # rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No Var")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'b', marker = "x", s = 3) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("No_Granulation/rm_and_residuals_{}.png".format(lines_no_cb[i]))
    plt.clf()

    # rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID line RVs")
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_rv_array)**2))/len(line_rv_array - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-40], -300, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals 
    axs[1].scatter(UTC_time, line_rv_array - GRASS_rv_array, color = 'b', marker = "x", s = 3) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("Granulation/LinebyLine/rm_and_residuals_cb_{}.png".format(lines[i]))
    plt.clf()

    # rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # residuals
    axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)    
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("Granulation/rm_and_residuals_cb_{}.png".format(lines[i]))
    plt.clf()