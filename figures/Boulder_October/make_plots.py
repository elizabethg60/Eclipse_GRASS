import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#read in data
#GRASS
grass_data = h5py.File("data/boulder_rv_regular.jld2", "r")
GRASS_rv  = grass_data["rv"][()]
rv_error_GRASS_cb  = grass_data["rv_error"][()]
grass_data_no_cb = h5py.File("data/boulder_rv_off.jld2", "r")
GRASS_no_cb  = grass_data_no_cb["rv"][()]
rv_error_GRASS_no_cb = grass_data_no_cb["rv_error"][()]
projected_no_cb = h5py.File("data/boulder_rv_projected.jld2", "r")
RV_list_no_cb  = projected_no_cb["RV_list_no_cb"][()]

#data 
data_time = np.loadtxt("data/Boulder_Data_bin.txt")[:, 0]
rv_obs = np.loadtxt("data/Boulder_Data_bin.txt")[:, 1] * (1.565*10**(-6))
time_stamps_data = []
time_julian = []
for i in data_time:
    dt = datetime.fromtimestamp(i) + timedelta(hours=4)
    time_stamps_data.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=39.995380 , longi=-105.262390 , alt=165.23, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs[47:-20])
rv_obs -= rv_obs[-1]

time_stamps_data = time_stamps_data[47:-20]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][47:-20]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def plot_line(time_stamps_data, rv_obs, model, model_label, GRASS, GRASS_label, LD, save, rv_error_GRASS_no_cb_array):
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(time_stamps_data, rv_obs, color = 'k', marker = "x", s = 18, label = "Real RVs")
    axs[0].plot(time_stamps_data, model, color = 'r', linewidth = 2, label = model_label)
    axs[0].plot(time_stamps_data, GRASS, color = 'b', linewidth = 2, label = GRASS_label)
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - model)**2))/len(rv_obs - model)),2)
    axs[0].text(time_stamps_data[-20], -300, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS)**2))/len(rv_obs - GRASS)),2)
    axs[0].text(time_stamps_data[-20], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    axs[0].text(time_stamps_data[-20], -700, "GRASS RV error {}".format(round(np.mean(rv_error_GRASS_no_cb_array),2)))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(time_stamps_data, rv_obs - model, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(time_stamps_data, rv_obs - GRASS, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("{}_{}.png".format(save, LD))
    plt.clf()

LD_arr = ["HD", "K300", "KSSD", "NIR"]
for i in range(0,len(LD_arr)):
    rv_error_GRASS_cb_array = grass_data[rv_error_GRASS_cb[i]][()][47:-20]
    rv_error_GRASS_no_cb_array = grass_data_no_cb[rv_error_GRASS_no_cb[i]][()][47:-20]
    
    GRASS_rv_array = jld2_read(grass_data, GRASS_rv, vb[47:-20], i)
    GRASS_no_cb_array = jld2_read(grass_data_no_cb, GRASS_no_cb, vb[47:-20], i)
    RV_list_no_cb_array = jld2_read(projected_no_cb, RV_list_no_cb, vb[47:-20], i)

    #rm curve 
    plot_line(time_stamps_data, rv_obs, RV_list_no_cb_array, "Weighted RVs - No CB", GRASS_no_cb_array, "Line RVs - No Var", LD_arr[i], "no_granulation", rv_error_GRASS_no_cb_array)
    # plot_line(time_stamps_data, rv_obs, 'RV_list_no_cb_array', "Weighted RVs - No CB", GRASS_rv_array, "Line RVs - Var", LD_arr[i], "granulation", rv_error_GRASS_cb_array)