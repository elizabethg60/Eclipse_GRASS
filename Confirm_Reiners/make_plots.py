import h5py
import math
import csv
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt 
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from astropy.time import Time

#data 
f=open("data/Reiners_Data.txt","r")
lines=f.readlines()[1:]
raw_rv = []
UTC_time = []
for x in lines:
    raw_rv.append(float(x.split()[4]))
for x in range(0,len(lines)):
    UTC_time.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f"))) #+ timedelta(seconds=float(50))
f.close()

#model
file = h5py.File("data/model_data_new.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
RV_list_cb_new  = file["RV_list_cb_new"][()]
intensity_list = file["intensity_list"][()]

#reiners model residuals and intensity 
with open('data/model2_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
with open('data/model1_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model1 = []
    residual_model1 = []
    for row in cf:
        time_residuals_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model1.append(float(row[1]))
with open('data/new_residuals_model1.csv') as f:
    cf = csv.reader(f)
    time_residuals_new1 = []
    residual_new1 = []
    for row in cf:
        time_residuals_new1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new1.append(float(row[1]))
with open('data/flux_data.csv') as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))    

raw_rv = np.array(raw_rv)
raw_rv -= raw_rv[-1]

RV_list_no_cb_array = np.array(RV_list_no_cb) 
RV_list_no_cb_array -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb)
RV_list_cb -= RV_list_cb[-1]

RV_list_cb_new = np.array(RV_list_cb_new)
RV_list_cb_new -= RV_list_cb_new[-1]

residual_slope = ((raw_rv[46] - RV_list_no_cb_array[46])-(raw_rv[0] - RV_list_no_cb_array[0]))/46
linear_correction = residual_slope*range(0,46) + (raw_rv[0] - RV_list_no_cb_array[0])

corrected_rv = list((raw_rv[0:46]) - linear_correction) + list(raw_rv[46:,])

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_no_cb_array, color = 'r', label = "Weighted RVs - No CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_no_cb_array)**2))/len(corrected_rv - RV_list_no_cb_array)),2)
# axs[0].text(UTC_time[-55], -400, "Weighted RVs - No CB RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals
# axs[1].scatter(UTC_time, corrected_rv - RV_list_no_cb_array, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model1, residual_model1, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_no_cb_mean.png")
# plt.show()

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Weighted RVs - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb)**2))/len(corrected_rv - RV_list_cb)),2)
# axs[0].text(UTC_time[-55], -400, "Weighted RVs - CB RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals       
# axs[1].scatter(UTC_time, corrected_rv - RV_list_cb, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_cb_mean.png")
# plt.show()

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb_new, color = 'r', label = "Weighted RVs - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb_new)**2))/len(corrected_rv - RV_list_cb_new)),2)
# axs[0].text(UTC_time[-55], -400, "Weighted RVs - CB RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals       
# axs[1].scatter(UTC_time, corrected_rv - RV_list_cb_new, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_new_cb_mean.png")
# plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(UTC_time, intensity_list/intensity_list[129], color = 'r', label = "model")
# ax1.scatter(time_flux_model1, flux_model1, color = 'k', label = "Reiners")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity")
# plt.legend()
# plt.savefig("intensity_mean.png")
# plt.show()

fig, axs = plt.subplots(2,2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]}, figsize=(10, 5))
axs[0,0].scatter(UTC_time, raw_rv, color = 'k', marker = "x", s = 15, label = "Reiners")
axs[0,0].plot(UTC_time, RV_list_no_cb_array, color = 'r', label = "My Model")
axs[0,0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0,0].set_xlabel("Time (UTC)")
axs[0,0].set_ylabel("RV [m/s]")
axs[0,0].legend()
#residuals
axs[1,0].scatter(UTC_time, raw_rv - RV_list_no_cb_array, color = 'r', marker = "x", s = 1)
axs[1,0].scatter(time_residuals_model1, residual_model1, color = 'k', marker = "x", s = 1)
axs[1,0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1,0].set_xlabel("Time (UTC)")
axs[1,0].set_ylabel("Residuals")
#intensity
axs[0,1].scatter(time_flux_model1, flux_model1, color = 'k', label = "Reiners", marker = "x", s = 15,)
axs[0,1].plot(UTC_time, intensity_list/intensity_list[129], color = 'r', label = "My Model")
axs[0,1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0,1].set_xlabel("Time (UTC)")
axs[0,1].set_ylabel("Relative Intensity")
# Remove the unused subplot in the right panel
fig.delaxes(axs[1, 1])
plt.savefig("2nd_yr_project.pdf")
plt.show()
