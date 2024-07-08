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
# #GRASS
# grass_data = h5py.File("data/300_neid_all_lines_rv_regular.jld2", "r")
# lines = grass_data["name"][()]
# GRASS_rv  = grass_data["rv"][()]
# rv_error_GRASS_cb  = grass_data["rv_error"][()]
# grass_data_no_cb = h5py.File("data/300_neid_all_lines_rv_off.jld2", "r")
# lines = grass_data_no_cb["name"][()]
# GRASS_no_cb  = grass_data_no_cb["rv"][()]
# rv_error_GRASS_no_cb = grass_data_no_cb["rv_error"][()]
#model
file = h5py.File("model_data_Kostogryz_LD_SSD.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]
RV_list_cb  = file["RV_list_cb"][()]
# intensity_list = file["intensity_list"][()]
#data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
# #line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_lines = line_data["name"][()]
# line_rv = line_data["rv"][()]
# rv_error_line = line_data["rv_error"][()]

rv_obs = np.array(rv_obs)[0:-25]
rv_obs -= rv_obs[-1]

UTC_time = UTC_time[0:-25]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

for i in range(0,len(line_lines)):

    # rv_error_GRASS_cb_array = grass_data[rv_error_GRASS_cb[i]][()][0:-25]
    # rv_error_GRASS_no_cb_array = grass_data_no_cb[rv_error_GRASS_no_cb[i]][()][0:-25]
    # rv_error_line_array = line_data[rv_error_line[i]][()][0:-25]
    
    # GRASS_rv_array = grass_data[GRASS_rv[i]][()][0:-25]
    # GRASS_rv_array = np.array(GRASS_rv_array + vb)
    # GRASS_rv_array -= GRASS_rv_array[-1]

    # line_rv_array = line_data[line_rv[i]][()][0:-25]
    # line_rv_array = np.array(line_rv_array + vb)
    # line_rv_array -= line_rv_array[-1]

    # GRASS_no_cb_array = grass_data_no_cb[GRASS_no_cb[i]][()][0:-25]
    # GRASS_no_cb_array = np.array(GRASS_no_cb_array + vb)
    # GRASS_no_cb_array -= GRASS_no_cb_array[-1]

    RV_list_no_cb_array = file[RV_list_no_cb[i]][()][0:-25]
    RV_list_no_cb_array = np.array(RV_list_no_cb_array + vb)
    RV_list_no_cb_array -= RV_list_no_cb_array[-1]

    RV_list_cb_array = file[RV_list_cb[i]][()][0:-25]
    RV_list_cb_array = np.array(RV_list_cb_array + vb)
    RV_list_cb_array -= RV_list_cb_array[-1]

    # #rm curve 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    # # axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID line RVs") 
    # axs[0].plot(UTC_time, RV_list_no_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - No CB")
    # # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("Time (UTC)", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # rms_model_no_cb = round(np.sqrt((np.nansum((line_rv_array - RV_list_no_cb_array)**2))/len(line_rv_array - RV_list_no_cb_array)),2)
    # axs[0].text(UTC_time[-40], 500, "Model RMS {}".format(rms_model_no_cb))
    # # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_no_cb_array)**2))/len(line_rv_array - GRASS_no_cb_array)),2)
    # # axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time, line_rv_array - RV_list_no_cb_array, color = 'r', marker = "x", s = 3) 
    # # axs[1].scatter(UTC_time, line_rv_array - GRASS_no_cb_array, color = 'b', marker = "x", s = 3)  
    # # axs[1].text(UTC_time[-70], 30, "GRASS - no CB avg error {}".format(round(np.mean(rv_error_GRASS_no_cb_array),2)))
    # axs[1].text(UTC_time[-70], 15, "Line RV avg error {}".format(round(np.mean(rv_error_line_array),2)))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array[120:-1])**2))/len(line_rv_array[120:-1])),2)
    # axs[1].text(UTC_time[-70], 0, "Out of Transit RMS {}".format(rms_grass_no_cb))
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("Time (UTC)", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("300_with_old_extinction/line_rm_and_residuals_{}.png".format(lines[i]))
    # plt.clf()

    #rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    axs[0].plot(UTC_time, RV_list_no_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - No CB")
    # axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No CB")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb_array)**2))/len(rv_obs - RV_list_no_cb_array)),2)
    axs[0].text(UTC_time[-40], 500, "Model RMS {}".format(rms_model_no_cb))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    # axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time, rv_obs - RV_list_no_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'b', marker = "x", s = 3) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("LD_SSD/no_cb/LD_SSD_rm_and_residuals_{}.png".format(line_lines[i]))
    #plt.show()
    plt.clf()

    # #rm curve w granulation 
    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    # axs[0].scatter(UTC_time, line_rv_array, color = 'y', marker = "x", s = 18, label = "NEID line RVs")
    # axs[0].plot(UTC_time, RV_list_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - CB")
    # # axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("Time (UTC)", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # # rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb_array)**2))/len(rv_obs - RV_list_cb_array)),2)
    # # axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    # # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    # # axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    # rms_model_no_cb = round(np.sqrt((np.nansum((line_rv_array - RV_list_cb_array)**2))/len(line_rv_array - RV_list_cb_array)),2)
    # axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    # # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS_rv_array)**2))/len(line_rv_array - GRASS_rv_array)),2)
    # # axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # # axs[1].scatter(UTC_time, rv_obs - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    # # axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    # axs[1].scatter(UTC_time, line_rv_array - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    # # axs[1].scatter(UTC_time, line_rv_array - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    # # axs[1].text(UTC_time[-80], 30, "GRASS - no CB avg error {}".format(round(np.mean(rv_error_GRASS_cb_array),2)))
    # axs[1].text(UTC_time[-80], 15, "Line RV avg error {}".format(round(np.mean(rv_error_line_array),2)))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array[120:-1])**2))/len(line_rv_array[120:-1])),2)
    # axs[1].text(UTC_time[-80], 0, "Out of Transit RMS {}".format(rms_grass_no_cb))
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("Time (UTC)", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("300_with_old_extinction/line_rm_and_residuals_cb_{}.png".format(lines[i]))
    # #plt.show()
    # plt.clf()

    #rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].plot(UTC_time, RV_list_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - CB")
    # axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb_array)**2))/len(rv_obs - RV_list_cb_array)),2)
    axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    # axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time, rv_obs - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)    
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("LD_SSD/cb/LD_SSD_rm_and_residuals_cb_{}.png".format(line_lines[i]))
    #plt.show()
    plt.clf()