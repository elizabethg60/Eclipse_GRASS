import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

GRASS_time = ["2023-10-14T15:26:45", "2023-10-14T15:28:07", "2023-10-14T15:29:30", "2023-10-14T15:30:53", "2023-10-14T15:32:15", "2023-10-14T15:33:38", "2023-10-14T15:35:01", "2023-10-14T15:36:23", "2023-10-14T15:37:46", "2023-10-14T15:39:09", "2023-10-14T15:40:31", "2023-10-14T15:41:54", "2023-10-14T15:43:17", "2023-10-14T15:44:39", "2023-10-14T15:46:02", "2023-10-14T15:47:25", "2023-10-14T15:48:47", "2023-10-14T15:50:10", "2023-10-14T15:51:33", "2023-10-14T15:52:56", "2023-10-14T15:54:18", "2023-10-14T15:55:41", "2023-10-14T15:57:04", "2023-10-14T15:58:26", "2023-10-14T15:59:49", "2023-10-14T16:01:12", "2023-10-14T16:02:34", "2023-10-14T16:03:57", "2023-10-14T16:05:20", "2023-10-14T16:06:42", "2023-10-14T16:08:05", "2023-10-14T16:09:28", "2023-10-14T16:10:50", "2023-10-14T16:12:13", "2023-10-14T16:13:36", "2023-10-14T16:14:58", "2023-10-14T16:16:21", "2023-10-14T16:17:44", "2023-10-14T16:19:06", "2023-10-14T16:20:29", "2023-10-14T16:21:52", "2023-10-14T16:23:15", "2023-10-14T16:24:37", "2023-10-14T16:26:00", "2023-10-14T16:27:23", "2023-10-14T16:28:45", "2023-10-14T16:30:08", "2023-10-14T16:31:31", "2023-10-14T16:32:53", "2023-10-14T16:34:16", "2023-10-14T16:35:39", "2023-10-14T16:37:01", "2023-10-14T16:38:24", "2023-10-14T16:39:47", "2023-10-14T16:41:09", "2023-10-14T16:42:32", "2023-10-14T16:43:55", "2023-10-14T16:45:17", "2023-10-14T16:46:40", "2023-10-14T16:48:03", "2023-10-14T16:49:25", "2023-10-14T16:50:48", "2023-10-14T16:52:11", "2023-10-14T16:53:33", "2023-10-14T16:54:56", "2023-10-14T16:56:19", "2023-10-14T16:57:42", "2023-10-14T16:59:04", "2023-10-14T17:00:27", "2023-10-14T17:01:50", "2023-10-14T17:03:12", "2023-10-14T17:04:35", "2023-10-14T17:05:58", "2023-10-14T17:07:20", "2023-10-14T17:08:43", "2023-10-14T17:10:06", "2023-10-14T17:11:28", "2023-10-14T17:12:51", "2023-10-14T17:14:14", "2023-10-14T17:15:36", "2023-10-14T17:16:59", "2023-10-14T17:18:22", "2023-10-14T17:19:44", "2023-10-14T17:21:07", "2023-10-14T17:22:30", "2023-10-14T17:23:52", "2023-10-14T17:25:15", "2023-10-14T17:26:38", "2023-10-14T17:28:01", "2023-10-14T17:29:23", "2023-10-14T17:30:46", "2023-10-14T17:32:09", "2023-10-14T17:33:31", "2023-10-14T17:34:54", "2023-10-14T17:36:17", "2023-10-14T17:37:39", "2023-10-14T17:39:02", "2023-10-14T17:40:25", "2023-10-14T17:41:47", "2023-10-14T17:43:10", "2023-10-14T17:44:33", "2023-10-14T17:45:55", "2023-10-14T17:47:18", "2023-10-14T17:48:41", "2023-10-14T17:50:03", "2023-10-14T17:51:26", "2023-10-14T17:52:49", "2023-10-14T17:54:11", "2023-10-14T17:55:34", "2023-10-14T17:56:57", "2023-10-14T17:58:20", "2023-10-14T17:59:42", "2023-10-14T18:01:05", "2023-10-14T18:02:28", "2023-10-14T18:03:50", "2023-10-14T18:05:13", "2023-10-14T18:06:36", "2023-10-14T18:07:58", "2023-10-14T18:09:21", "2023-10-14T18:10:44", "2023-10-14T18:12:06", "2023-10-14T18:13:29", "2023-10-14T18:14:52", "2023-10-14T18:16:14", "2023-10-14T18:17:37", "2023-10-14T18:19:00", "2023-10-14T18:20:22", "2023-10-14T18:21:45", "2023-10-14T18:23:08", "2023-10-14T18:24:30", "2023-10-14T18:25:53", "2023-10-14T18:27:16", "2023-10-14T18:28:38", "2023-10-14T18:30:01", "2023-10-14T18:31:24", "2023-10-14T18:32:47", "2023-10-14T18:34:09", "2023-10-14T18:35:32", "2023-10-14T18:36:55", "2023-10-14T18:38:17", "2023-10-14T18:39:40", "2023-10-14T18:41:03", "2023-10-14T18:42:25", "2023-10-14T18:43:48", "2023-10-14T18:45:11", "2023-10-14T18:46:33", "2023-10-14T18:47:56", "2023-10-14T18:49:19", "2023-10-14T18:50:41", "2023-10-14T18:52:04", "2023-10-14T18:53:27", "2023-10-14T18:54:49", "2023-10-14T18:56:12", "2023-10-14T18:57:35", "2023-10-14T18:58:57", "2023-10-14T19:00:20", "2023-10-14T19:01:43", "2023-10-14T19:03:06"]

#read in data
#GRASS
grass_data = h5py.File("neid_all_lines_rv.jld2", "r")
lines = grass_data["name"][()]
grass_rv  = grass_data["rv"][()]
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]
RV_list_cb  = file["RV_list_cb"][()]
intensity_list = file["intensity_list"][()]
#data 
data = pd.read_csv("NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
GRASS_time_st = []
for i in GRASS_time:
    GRASS_time_st.append(datetime.strptime(i, "%Y-%m-%dT%H:%M:%S"))
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs)
rv_obs -= rv_obs[-1]

RV_list_no_cb = np.array(RV_list_no_cb + vb)
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb + vb)
RV_list_cb -= RV_list_cb[-1]

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', linewidth = 2, label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
axs[0].legend(fontsize=12)
rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb)**2))/len(rv_obs - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(UTC_time, (rv_obs) - RV_list_no_cb, color = 'r', marker = "x", s = 3) 
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)", fontsize=12)
axs[1].set_ylabel("Residuals", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("rm_and_residuals.png")
plt.show()

for i in range(1,len(lines)):
    GRASS_rv = grass_data[grass_rv[i]][()]
    GRASS_rv = np.array(GRASS_rv + vb)
    GRASS_rv -= GRASS_rv[-1]

    #rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', linewidth = 2, label = "Model - No CB")
    axs[0].plot(GRASS_time_st, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb)**2))/len(rv_obs - RV_list_no_cb)),2)
    axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv)**2))/len(rv_obs - GRASS_rv)),2)
    axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time, (rv_obs) - RV_list_no_cb, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, (GRASS_rv) - RV_list_no_cb, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("GRASS/rm_and_residuals_{}.png".format(lines[i]))
    plt.show()

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
axs[0].plot(UTC_time, RV_list_cb, color = 'r', linewidth = 3, label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
rms_model_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb)**2))/len(rv_obs - RV_list_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_cb))
axs[0].legend(fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(UTC_time, (rv_obs) - RV_list_cb, color = 'r', marker = "x", s = 3)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)", fontsize=12)
axs[1].set_ylabel("Residuals", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("rm_and_residuals_cb.png")
plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(UTC_time, intensity_list, label = "Model") 
# #ax1.scatter(UTC_time, i_test, color = 'g', label = "GRASS")  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# #plt.legend()
# plt.savefig("intensity.png")
# plt.show()