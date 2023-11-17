import h5py
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.time import Time
import numpy as np
from barycorrpy import get_BC_vel, exposure_meter_BC_vel
import jdcal
mpl = plt.matplotlib 

#rm curve
f = h5py.File("rv_intensity.jld2", "r")
RV_list = f["RV_list"][()]
intensity_list = f["intensity_list"][()]
date_strings = f["timestamps"][()]
fig = plt.figure()
ax1 = fig.add_subplot()
time_stamps = []
for i in date_strings:
    time_stamps.append(datetime.strptime(i, "%Y-%m-%dT%H:%M:%S"))
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
vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=51.54548, longi=9.905548, alt=150, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
ax1.scatter(time, raw_rv + vb, label = "Reiners corrected RVs")
ax1.scatter(time_stamps, RV_list+ (raw_rv[-1] + vb[-1] - RV_list[-1]), label = "Model")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("RV [m/s]")
plt.legend()
plt.savefig("rm_curve.png")
plt.show()
#residuals
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(time_stamps, (raw_rv + vb) - (RV_list+ (raw_rv[-1] + vb[-1] - RV_list[-1])))  
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("residuals") 
plt.savefig("residuals.png")
plt.show()

#intensity
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(time_stamps, intensity_list)  
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("Time (UTC)")
ax1.set_ylabel("Relative Intensity") 
plt.savefig("intensity.png")
plt.show()

#projected solar velocities at each timestamp for eclipse movie
for i in range(1,160):
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

