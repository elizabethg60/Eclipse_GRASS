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
grass_data = h5py.File("neid_all_lines_rv_LD.jld2", "r")
lines = grass_data["name"][()]
grass_rv  = grass_data["rv"][()]
#model
file = h5py.File("model_data.jld2", "r")
RV_list_model = file["RV_list_no_cb"][()]

#data 
data = pd.read_csv("neid_april_data.csv")
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
    GRASS_rv = grass_data[grass_rv[i]][()][24:,]
    GRASS_rv = np.array(GRASS_rv + vb)
    GRASS_rv -= GRASS_rv[-1]

    RV_list_no_cb = file[RV_list_model[i]][()][24:,]
    RV_list_no_cb = np.array(RV_list_no_cb + vb)
    RV_list_no_cb -= RV_list_no_cb[-1]

    #rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', linewidth = 2, label = "Model - No CB")
    axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb)**2))/len(rv_obs - RV_list_no_cb)),2)
    axs[0].text(UTC_time[-60], -300, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv)**2))/len(rv_obs - GRASS_rv)),2)
    axs[0].text(UTC_time[-60], -400, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
 
    #residuals
    axs[1].scatter(UTC_time, rv_obs - RV_list_no_cb, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, rv_obs - GRASS_rv, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("GRASS/rm_and_residuals_{}.png".format(lines[i]))
    #plt.show()
    plt.clf()

    #clouds time index 2,4,11,15