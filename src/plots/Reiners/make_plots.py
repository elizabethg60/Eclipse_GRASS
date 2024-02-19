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

#read in data
#model
#file = h5py.File("Resolution/model_data100_20sub.jld2", "r")
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
RV_list_cb_new  = file["RV_list_cb_new"][()]
intensity_list = file["intensity_list"][()]
RA_list = file["RA_list"][()]
dec_list = file["dec_list"][()]
vel_no_cb = file["vel_no_cb"][()]
vel_cb = file["vel_cb"][()]

#data 
f=open("Reiners_Data.txt","r")
lines=f.readlines()[1:]
time=[]
raw_rv = []
time_julian = []
UTC_time = []
for x in lines:
    dt = datetime.strptime("2015-03-20 {}".format((x.split()[1])), "%Y-%m-%d %H:%M:%S.%f")
    time.append(dt)
    time_julian.append((Time(dt)).jd)
    raw_rv.append(float(x.split()[4]))
for x in range(0,len(lines)):
    UTC_time.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f") + timedelta(seconds=float(50))))
f.close()
#reiners model one residuals and intensity 
with open('model1_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model1 = []
    residual_model1 = []
    for row in cf:
        time_residuals_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model1.append(float(row[1]))
with open('model2_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
with open('new_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_new = []
    residual_new = []
    for row in cf:
        time_residuals_new.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new.append(float(row[1]))
with open('flux_data.csv') as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))     

raw_rv = np.array(raw_rv)
raw_rv -= raw_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb) 
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb)
RV_list_cb -= RV_list_cb[-1]

RV_list_cb_new = np.array(RV_list_cb_new)
RV_list_cb_new -= RV_list_cb_new[-1]

residual_slope = ((raw_rv[46] - RV_list_cb[46])-(raw_rv[0] - RV_list_cb[0]))/46
linear_correction = residual_slope*range(0,46) + (raw_rv[0] - RV_list_cb[0])

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_no_cb[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_no_cb[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model1, residual_model1, color = 'k', marker = "x", s = 1)
axs[1].scatter(time_residuals_new, residual_new, color = 'b', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals       
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb_new, color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb_new[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb_new[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_cb_new.png")
plt.show()

#intensity
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(UTC_time, intensity_list/intensity_list[129], color = 'r', label = "model")
ax1.scatter(time_flux_model1, flux_model1, color = 'k', label = "Reiners")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("Relative Intensity")
plt.legend()
plt.savefig("intensity.png")
plt.show()

# #projected solar velocities at each timestamp for eclipse movie
# for i in range(1,len(vel_no_cb)):
#     vel = file[vel_no_cb[i]][()]
#     ra = file[RA_list[i]][()]
#     dec = file[dec_list[i]][()]

#     cnorm = mpl.colors.Normalize(np.min(vel), np.max(vel))
#     colors = mpl.cm.seismic(cnorm(vel))
#     pcm = plt.pcolormesh(ra, dec, vel, cmap="seismic",vmin=-2000, vmax=2000)
#     cb = plt.colorbar(pcm, norm=cnorm, ax=plt.gca())
#     plt.gca().invert_xaxis()
#     cb.set_label("projected velocity (m/s)")
#     plt.xlabel("RA")
#     plt.ylabel("dec")
#     plt.savefig("movie/projected_vel_{}.png".format(i))
#     plt.clf()

# #barycentric correction test: Sun center vs surface
# vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=51.560583 , longi=9.944333, alt=201, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, vb, color = 'g', label = "BC using center")
# plt.scatter(time, bc_list, color = 'r', label = "BC using surface")
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Barycentric Correction")
# plt.legend()
# plt.savefig("barycentric_correction.png")
# plt.show()

# #extiction test: differential vs absolute 
# f=open("Reiners_Data.txt","r")
# lines=f.readlines()[1:]
# airmass = []
# for x in range(0,len(lines)):
#    airmass.append(float((lines[x].split()[2])))
# f.close()
# ext_coefficient = list(np.linspace(0.4, 0.1, num = len(airmass)))
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, airmass, color = 'b', label = "absolute airmass")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Airmass")
# diff_airmass1 = []
# diff_airmass2 = []
# for i in range(0,len(airmass_matrix)):
#     diff_airmass1.append(1 / file[airmass_matrix[i]][()][1][1])
#     diff_airmass2.append(1 / file[airmass_matrix[i]][()][30][10])
# plt.scatter(time, diff_airmass1, color = 'r', label = "airmass for random patch")
# plt.scatter(time, diff_airmass2, color = 'g', label = "airmass for random patch")
# plt.legend()
# plt.savefig("airmass_test")
# plt.show()
