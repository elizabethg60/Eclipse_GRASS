import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel
import matplotlib.lines as mlines
import matplotlib.colors as mcolors

# projected_SSD_4parameter = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB/data/projected_SSD_4parameter_gpu_ext_CB_optim.jld2", "r")
projected_SSD_4parameter = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/Optim/CB_MF/data/projected_SSD_4parameter_gpu_ext_CB_MF_optim.jld2", "r")
projected_RV_SSD_4parameter = projected_SSD_4parameter["RV_list_no_cb"][()]

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_lines = line_data["name"][()]
line_rv  = line_data["rv"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
UTC_time = UTC_time[0:-28]

wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
norm = mcolors.Normalize(vmin=np.min(wavelength), vmax=np.max(wavelength))
cmap = plt.get_cmap('coolwarm')

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + vb)
    array -= array[-1]
    return array

def bin_array(arr, bin_size):
    arr = np.array(arr)
    n_bins = len(arr) // bin_size
    return np.array([np.nanmean(arr[i*bin_size:(i+1)*bin_size]) for i in range(n_bins)])

def bin_array_time(arr, bin_size):
    arr = np.array(arr)
    n_bins = len(arr) // bin_size
    return np.array([(arr[i*bin_size:(i+1)*bin_size])[0] for i in range(n_bins)])

line_list = ["FeI_5250.2","FeI_5250.5","FeI_5379","TiII_5381","FeI_5382","FeI_5383","MnI_5432","FeI_5432","FeI_5434","NiI_5435","FeI_5436.3","FeI_5436.6","FeI_5576","NiI_5578","FeII_6149","FeI_6151","CaI_6169.0","CaI_6169.5","FeI_6170","FeI_6173","FeI_6301","FeI_6302"]

fig, axs = plt.subplots(1, figsize=(8, 4))

res_arr = []
for i, value in enumerate(wavelength):
    color = cmap(norm(value))     
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    projected_RV_SSD_4parameter_array = jld2_read(projected_SSD_4parameter, projected_RV_SSD_4parameter, vb, i)
    axs.scatter(UTC_time, line_rv_array - projected_RV_SSD_4parameter_array, color = color, marker = "x", s = 8) 
    res_arr.append(line_rv_array - projected_RV_SSD_4parameter_array)
# means = [sum(col) / len(col) for col in zip(*res_arr)]
# axs.scatter(UTC_time, means, color = "k", marker = "x", s = 8) 
axs.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs.set_xlabel("Time (UTC)", fontsize=15)
axs.set_ylabel("Residuals (m/s)", fontsize=15) 
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   
cbar = plt.colorbar(sm, ax=axs)
cbar.set_label("Wavelength (Å)", fontsize=14)
cbar.ax.tick_params(labelsize=11)
plt.tick_params(axis='x', which='both', top=True, labeltop=False, labelsize = 13)
plt.tick_params(axis='y', labelsize = 13)
plt.savefig("residuals_MF.pdf", bbox_inches="tight")
plt.clf()

  