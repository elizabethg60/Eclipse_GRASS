import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#GRASS_rv: RV calculated from line profiles + granulation affects on
#GRASS_no_cb: RV calculated from line profiles + granulation affects off
#RV_list_no_cb: RV calculated form weighted projected velocities + no extra CB
#RV_list_cb: RV calculated form weighted projected velocities + with extra CB

#read in data
#GRASS
grass_data = h5py.File("neid_all_lines_rv_regular.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
grass_data_no_cb = h5py.File("neid_all_lines_rv_off.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]
RV_list_cb  = file["RV_list_cb"][()]

#data 
data = pd.read_csv("data/neid_april_data.csv")
rv_obs = list(data["ccfrvmod"][124:-10]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][124:-10]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs)
rv_obs -= rv_obs[-1]

for i in range(1,len(lines)):
    GRASS_rv_array = grass_data[GRASS_rv[i]][()][24:,]
    GRASS_rv_array = np.array(GRASS_rv_array + vb)
    GRASS_rv_array -= GRASS_rv_array[-1]

    GRASS_no_cb_array = grass_data_no_cb[GRASS_no_cb[i]][()][24:,]
    GRASS_no_cb_array = np.array(GRASS_no_cb_array + vb)
    GRASS_no_cb_array -= GRASS_no_cb_array[-1]

    RV_list_no_cb_array = file[RV_list_no_cb[i]][()][24:,]
    RV_list_no_cb_array = np.array(RV_list_no_cb_array + vb)
    RV_list_no_cb_array -= RV_list_no_cb_array[-1]

    RV_list_cb_array = file[RV_list_cb[i]][()][24:,]
    RV_list_cb_array = np.array(RV_list_cb_array + vb)
    RV_list_cb_array -= RV_list_cb_array[-1]

    #rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, RV_list_no_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - No CB")
    axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No CB")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb_array)**2))/len(rv_obs - RV_list_no_cb_array)),2)
    axs[0].text(UTC_time[-60], -200, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    axs[0].text(UTC_time[-60], -300, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time, rv_obs - RV_list_no_cb_array, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("No_Granulation/rm_and_residuals_{}.png".format(lines_no_cb[i]))
    #plt.show()
    plt.clf()

    #rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, RV_list_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - CB")
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb_array)**2))/len(rv_obs - RV_list_cb_array)),2)
    axs[0].text(UTC_time[-60], -200, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-60], -300, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time, rv_obs - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("Granulation/rm_and_residuals_cb_{}.png".format(lines[i]))
    #plt.show()
    plt.clf()

    #clouds time index 2,4,11,15