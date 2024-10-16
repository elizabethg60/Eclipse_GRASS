import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#Boulder
projected_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/Boulder_October/data/boulder_rv_projected.jld2", "r")
RV_list_no_cb  = projected_no_cb["RV_list_no_cb"][()]
#data 
data_time = np.loadtxt("/storage/home/efg5335/work/Eclipse_GRASS/figures/Boulder_October/data/Boulder_Data_bin.txt")[:, 0]
rv_obs_LHR = np.loadtxt("/storage/home/efg5335/work/Eclipse_GRASS/figures/Boulder_October/data/Boulder_Data_bin.txt")[:, 1] * (1.565*10**(-6))
time_stamps_data = []
time_julian = []
for i in data_time:
    dt = datetime.fromtimestamp(i) + timedelta(hours=4)
    time_stamps_data.append(dt)
    time_julian.append((Time(dt)).jd)
vb_LHR, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=39.995380 , longi=-105.262390 , alt=165.23, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
rv_obs_LHR = np.array(rv_obs_LHR[47:-20])
rv_obs_LHR -= rv_obs_LHR[-1]
time_stamps_data = time_stamps_data[47:-20]

#NEID
grass_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KHD/data/neid_all_lines_rv_regular_KHD.jld2", "r")
GRASS_rv  = grass_data["rv"][()]
file = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KHD/data/neid_october_N_50_KHD.jld2", "r")
RV_list_no_cb_neid = file["RV_list_no_cb"][()]
#line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_rv  = line_data["rv"][()]
#data
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
vb_neid, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
UTC_time = UTC_time[0:-25]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][47:-20]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def jld2_read_neid(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    array -= array[-1]
    return array

RV_list_no_cb_array = jld2_read(projected_no_cb, RV_list_no_cb, vb_LHR[47:-20], 0)
GRASS_rv_array = jld2_read_neid(grass_data, GRASS_rv, vb_neid, 8)
line_rv_array = jld2_read_neid(line_data, line_rv, vb_neid, 8)
RV_list_no_cb_array_neid = jld2_read_neid(file, RV_list_no_cb_neid, vb_neid, 8)

# Create a figure and subplots
fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [10, 4]}, figsize=(10, 5))

# Plot the first figure
axs[0, 0].scatter(UTC_time, line_rv_array, color = 'k', marker = "x", s = 18, label = "NEID RVs")
axs[0, 0].set_title('NEID - 543.45 nm')
axs[0, 0].plot(UTC_time, RV_list_no_cb_array_neid, color = 'r', linewidth = 2, label = "Projected RVs")
axs[0, 0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
axs[0, 0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0, 0].set_xlabel("Time (UTC)", fontsize=12)
axs[0, 0].set_ylabel("RV [m/s]", fontsize=12)
axs[0, 0].legend(fontsize=12)
# Plot the residuals for the first figure
axs[1, 0].scatter(UTC_time, line_rv_array - GRASS_rv_array, color = 'b', marker = "x", s = 3)
axs[1, 0].scatter(UTC_time, line_rv_array - RV_list_no_cb_array_neid, color = 'r', marker = "x", s = 3)
axs[1, 0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1, 0].set_xlabel("Time (UTC)", fontsize=12)
axs[1, 0].set_ylabel("Residuals", fontsize=12) 
# Plot the second figure
axs[0, 1].scatter(time_stamps_data, rv_obs_LHR, color = 'k', marker = "x", s = 18, label = "LHR RVs")
axs[0, 1].plot(time_stamps_data, RV_list_no_cb_array, color = 'r', linewidth = 2, label = "Projected RVs")
axs[0, 1].set_title('LHR - 1565.2 nm')
axs[0, 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0, 1].legend(fontsize=12)
# Plot the residuals for the second figure
axs[1, 1].scatter(time_stamps_data, rv_obs_LHR - RV_list_no_cb_array, color = 'r', marker = "x", s = 3)
axs[1, 1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1, 1].set_xlabel("Time (UTC)", fontsize=12)
# Adjust the layout
plt.tight_layout()
plt.savefig("RM_curve_comp.png")
