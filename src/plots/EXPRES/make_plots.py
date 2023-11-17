import h5py
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
mpl = plt.matplotlib 

#rm curve
f = h5py.File("rv_intensity.jld2", "r")
RV_list = f["RV_list"][()]
intensity_list = f["intensity_list"][()]
date_strings = f["timestamps"][()]
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
time_stamps = []
for i in date_strings:
    time_stamps.append(datetime.strptime(i, "%Y-%m-%dT%H:%M:%S"))
axs[0].scatter(time_stamps[20:-100], RV_list[20:-100], color = 'k', label = "Model")
data = pd.read_csv("EXPRES_data.csv")
time_stamps_data = []
for i in data["tobs"]:
    time_stamps_data.append(datetime.strptime(i, "%Y-%m-%d %H:%M:%S"))
offset = list(data["rv"])[-1] - RV_list[-1]
axs[0].scatter(time_stamps_data[20:-100], data["rv"][20:-100] - offset, color = 'r', marker = "x", s = 5, label = "EXPRESS RVs") 
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(time_stamps[20:-100], (data["rv"][20:-100]-offset) - RV_list[20:-100], color = 'k', marker = "x", s = 1)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals") 
plt.savefig("rm_and_residuals.png")
plt.show()

#intensity
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(time_stamps, intensity_list)  
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("Relative Intensity") 
plt.savefig("intensity.png")
plt.show()

#projected solar velocities at each timestamp for eclipse movie
for i in range(1,len(time_stamps)):
    f = h5py.File("data/timestamp_{}.jld2".format(i), "r")
    projected_velocities = f["projected_velocities"][()]
    ra = f["ra"][()]
    dec = f["dec"][()]
    ra_moon = f["ra_moon"][()]
    dec_moon = f["dec_moon"][()]
    mu_grid = f["mu_grid"][()]
    LD_all = f["LD_all"][()]
    date_strings = f["timestamp"][()]

    #projected velocity rotation 
    cnorm = mpl.colors.Normalize(np.min(projected_velocities), np.max(projected_velocities))
    colors = mpl.cm.seismic(cnorm(projected_velocities))
    pcm = plt.pcolormesh(ra, dec, projected_velocities, cmap="seismic",vmin=-2000, vmax=2000)
    cb = plt.colorbar(pcm, norm=cnorm, ax=plt.gca())
    plt.gca().invert_xaxis()
    cb.set_label("projected velocity (m/s)")
    plt.xlabel("RA")
    plt.ylabel("dec")
    plt.savefig("movie/projected_vel_{}.png".format(i))
    plt.clf()