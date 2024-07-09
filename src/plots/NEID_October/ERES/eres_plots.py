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

#read in data
#GRASS
grass_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/neid_all_lines_rv_regular_ext.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
grass_data_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/neid_all_lines_rv_off_ext.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
#model
file = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/model_data_ext.jld2", "r")
intensity_list = file["intensity_list"][()]
#data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
#line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/neid_RVlinebyline_ext.jld2", "r")
line_lines = line_data["name"][()]
line_rv  = line_data["rv"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs)[0:-25]
rv_obs -= rv_obs[-1]

UTC_time = UTC_time[0:-25]

# for i in [5]:
#     GRASS_rv_array = grass_data[GRASS_rv[i]][()][0:-25]
#     GRASS_rv_array = np.array(GRASS_rv_array + vb)
#     GRASS_rv_array -= GRASS_rv_array[-1]

#     line_rv_array = line_data[line_rv[i]][()][0:-25]
#     line_rv_array = np.array(line_rv_array + vb)
#     line_rv_array -= line_rv_array[-1]

#     GRASS_no_cb_array = grass_data_no_cb[GRASS_no_cb[i]][()][0:-25]
#     GRASS_no_cb_array = np.array(GRASS_no_cb_array + vb)
#     GRASS_no_cb_array -= GRASS_no_cb_array[-1]

#     #rm curve 
#     fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
#     axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
#     axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
#     axs[0].set_ylabel("RV [m/s]", fontsize=12)
#     axs[0].legend(fontsize=12)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     #residuals
#     axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[1].set_xlabel("10/14/2023 Time (UTC)  ", fontsize=12)
#     axs[1].set_ylabel("Residuals", fontsize=12) 
#     axs[1].get_yaxis().set_ticks([])
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     plt.savefig("plot_one.png")
#     plt.clf()

    # #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'r', linewidth = 2, label = "GRASS - No CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    # print("image two")
    # print(rms_grass_no_cb)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("plot_two.png")
    # plt.clf()

    # #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'r', linewidth = 2, label = "GRASS - No CB")
    # axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS - CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    # print("image three")
    # print(rms_grass_no_cb)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("plot_three.png")
    # plt.clf()

    # #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID Line RVs")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("10/14/2023 Time (UTC)  ", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # axs[1].get_yaxis().set_ticks([])
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("plot_four.png")
    # plt.clf()

    #     #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID Line RVs")
    # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'r', linewidth = 2, label = "GRASS - No CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_no_cb_array)**2))/len(line_rv_array - GRASS_no_cb_array)),2)
    # print("image five")
    # print(rms_grass_no_cb)
    # axs[0].legend(fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time, line_rv_array - GRASS_no_cb_array, color = 'r', marker = "x", s = 3)  
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("10/14/2023 Time (UTC)  ", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # axs[1].get_yaxis().set_ticks([])
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("plot_five.png")
    # plt.clf()

    # #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID Line RVs")
    # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'r', linewidth = 2, label = "GRASS - No CB")
    # axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS - CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("10/14/2023 Time (UTC) ", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_rv_array)**2))/len(line_rv_array - GRASS_rv_array)),2)
    # print("image six")
    # print(rms_grass_no_cb)
    # axs[0].legend(fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time, line_rv_array - GRASS_no_cb_array, color = 'r', marker = "x", s = 3)  
    # axs[1].scatter(UTC_time, line_rv_array - GRASS_rv_array, color = 'b', marker = "x", s = 3) 
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("10/14/2023 Time (UTC)  ", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # axs[1].get_yaxis().set_ticks([])
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("plot_six.png")
    # plt.clf()

RMS_no_CB = []
RMS_CB = []
for i in range(0,len(lines)):
    GRASS_rv_array = grass_data[GRASS_rv[i]][()][0:-25]
    GRASS_rv_array = np.array(GRASS_rv_array + vb)
    GRASS_rv_array -= GRASS_rv_array[-1]

    line_rv_array = line_data[line_rv[i]][()][0:-25]
    line_rv_array = np.array(line_rv_array + vb)
    line_rv_array -= line_rv_array[-1]

    GRASS_no_cb_array = grass_data_no_cb[GRASS_no_cb[i]][()][0:-25]
    GRASS_no_cb_array = np.array(GRASS_no_cb_array + vb)
    GRASS_no_cb_array -= GRASS_no_cb_array[-1]

    RMS_no_CB.append(round(np.sqrt((np.nansum((line_rv_array - GRASS_no_cb_array)**2))/len(line_rv_array - GRASS_no_cb_array)),2))
    RMS_CB.append(round(np.sqrt((np.nansum((line_rv_array - GRASS_rv_array)**2))/len(line_rv_array - GRASS_rv_array)),2))

plt.figure(figsize=(12, 6))
plt.scatter(line_lines, RMS_no_CB, color = 'r', label = 'No CB')
plt.scatter(line_lines, RMS_CB, color = 'b', label = 'CB')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("plot_seven.png", bbox_inches='tight')