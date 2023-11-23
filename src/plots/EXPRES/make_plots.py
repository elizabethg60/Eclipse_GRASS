import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta

#read in data
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()][20:-100]
RV_list_cb  = file["RV_list_cb"][()][20:-100]
intensity_list = file["intensity_list"][()][20:-100]
RA_list = file["RA_list"][()]
dec_list = file["dec_list"][()]
vel_no_cb = file["vel_no_cb"][()]
vel_cb = file["vel_cb"][()]
#data 
data = pd.read_csv("EXPRES_data.csv")
rv_obs = list(data["rv"][20:-100])
UTC_time = []
for i in data["tobs"][20:-100]:
    UTC_time.append(datetime.strptime(i, "%Y-%m-%d %H:%M:%S"))

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, rv_obs, color = 'r', marker = "x", s = 5, label = "EXPRES RVs") 
offset = rv_obs[-1] - RV_list_no_cb[-1]
axs[0].scatter(UTC_time, RV_list_no_cb + offset, color = 'k', s = 5, label = "Model - No CB")
offset_cb = rv_obs[-1] - RV_list_cb[-1]
axs[0].scatter(UTC_time, RV_list_cb + offset_cb, color = 'grey', s = 5, label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time, (rv_obs-offset) - RV_list_no_cb, color = 'k', marker = "x", s = 1)  
axs[1].scatter(UTC_time, (rv_obs-offset_cb) - RV_list_cb, color = 'grey', marker = "x", s = 1)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals") 
plt.savefig("rm_and_residuals.png")
plt.show()

#intensity
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(UTC_time, intensity_list)  
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("Relative Intensity") 
plt.savefig("intensity.png")
plt.show()

#projected solar velocities at each timestamp for eclipse movie
for i in range(1,len(vel_no_cb)):
    vel = file[vel_no_cb[i]][()]
    ra = file[RA_list[i]][()]
    dec = file[dec_list[i]][()]

    cnorm = mpl.colors.Normalize(np.min(vel), np.max(vel))
    colors = mpl.cm.seismic(cnorm(vel))
    pcm = plt.pcolormesh(ra, dec, vel, cmap="seismic",vmin=-2000, vmax=2000)
    cb = plt.colorbar(pcm, norm=cnorm, ax=plt.gca())
    plt.gca().invert_xaxis()
    cb.set_label("projected velocity (m/s)")
    plt.xlabel("RA")
    plt.ylabel("dec")
    plt.savefig("movie/projected_vel_{}.png".format(i))
    plt.clf()