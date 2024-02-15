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
data_time = np.loadtxt("Boulder_Data_bin.txt")[:, 0]
data_rv = np.loadtxt("Boulder_Data_bin.txt")[:, 1] * (1.565*10**(-6))
time_stamps_data = []
time_julian = []
for i in data_time:
    dt = datetime.fromtimestamp(i) + timedelta(hours=4)
    time_stamps_data.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=39.995380 , longi=-105.262390 , alt=165.23, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

data_rv = np.array(data_rv)
data_rv -= data_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb + vb)
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb + vb)
RV_list_cb -= RV_list_cb[-1]

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(time_stamps_data[47:-20], data_rv[47:-20], color = 'k', marker = "x", s = 18, label = "Heterodyne RVs") 
axs[0].plot(time_stamps_data[47:-20], RV_list_no_cb[47:-20], color = 'r', linewidth = 3, label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
axs[0].legend(fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(time_stamps_data[47:-20], (data_rv[47:-20]) - RV_list_no_cb[47:-20], color = 'k', marker = "x", s = 3) 
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)", fontsize=12)
axs[1].set_ylabel("Residuals", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(time_stamps_data[47:-20], data_rv[47:-20], color = 'k', marker = "x", s = 18, label = "Heterodyne RVs") 
axs[0].plot(time_stamps_data[47:-20], RV_list_cb[47:-20], color = 'r', linewidth = 3, label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
axs[0].legend(fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(time_stamps_data[47:-20], (data_rv[47:-20]) - RV_list_cb[47:-20], color = 'k', marker = "x", s = 3)  
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
# ax1.scatter(time_stamps_data, intensity_list)  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# plt.savefig("intensity.png")
# plt.show()

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