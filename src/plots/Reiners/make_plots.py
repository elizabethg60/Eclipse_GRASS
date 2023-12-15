import h5py
import jdcal
import csv
import numpy as np
import pandas as pd
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt 
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#read in data
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
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
for x in lines:
    dt = datetime.strptime("2015-03-20 {}".format((x.split()[1])), "%Y-%m-%d %H:%M:%S.%f")
    time.append(dt)
    time_julian.append((Time(dt)).jd)
    raw_rv.append(float(x.split()[4]))
f.close()
#reiners model one residuals and intensity 
with open('model1_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model1 = []
    residual_model1 = []
    for row in cf:
        time_residuals_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model1.append(float(row[1]))
with open('flux_data.csv') as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))       

#rm curve 
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=51.560583 , longi=9.944333, alt=201, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
# v_obs = (raw_rv + vb) - (raw_rv + vb)[-1]
# axs[0].scatter(time, v_obs, color = 'r', marker = "x", s = 5, label = "Reiners corrected RVs")
# offset1 = (v_obs[-1] - RV_list_no_cb[-1])
# axs[0].scatter(time, RV_list_no_cb + offset1, color = 'k', s = 5, label = "Model - No CB")
# offset2 = (v_obs[-1] - RV_list_cb[-1])
# axs[0].scatter(time, RV_list_cb + offset2, color = 'grey', s = 5, label = "Model - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals
# axs[1].scatter(time, (v_obs-offset1) - RV_list_no_cb, color = 'k', marker = "x", s = 1)  
# axs[1].scatter(time, (v_obs-offset2) - RV_list_cb, color = 'grey', marker = "x", s = 1)  
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals") 
# plt.savefig("rm_and_residuals_bc.png")
# plt.show()

fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(time, raw_rv, color = 'k', marker = "x", s = 5, label = "Reiners RVs")
offset1 = (raw_rv[-1] - (RV_list_no_cb - vb)[-1])
axs[0].plot(time, (RV_list_no_cb - vb + offset1), color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(time, (raw_rv-offset1) - (RV_list_no_cb - vb), color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model1, residual_model1, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals") 
plt.savefig("rm_and_residuals.png")
plt.show()

#intensity
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(time, intensity_list/intensity_list[129]) 
ax1.scatter(time_flux_model1, flux_model1)  
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("Relative Intensity") 
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