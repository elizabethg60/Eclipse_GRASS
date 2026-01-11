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

SSD_quadratic = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD_4Parameter/LD/data/neid_all_lines_rv_on_SSD_4parameter_gpu.jld2", "r")
rv_gaussian = SSD_quadratic["rv"][()]
rv_symmetric_lsf = SSD_quadratic["rv_symmetric_lsf"][()]
rv_asymmetric_lsf = SSD_quadratic["rv_asymmetric_lsf"][()]

# data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

# line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r")
line_lines = line_data["name"][()]
line_rv  = line_data["rv"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-28], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
UTC_time = UTC_time[0:-28]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + vb)
    array -= array[-1]
    return array

line_list = ["FeI_5250.2","FeI_5250.5","FeI_5379","TiII_5381","FeI_5382","FeI_5383","MnI_5432","FeI_5432","FeI_5434","NiI_5435","FeI_5436.3","FeI_5436.6","FeI_5576","NiI_5578","FeII_6149","FeI_6151","CaI_6169.0","CaI_6169.5","FeI_6170","FeI_6173","FeI_6301","FeI_6302"]

for i in range(0,len(line_list)):    
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    rv_gaussian_array = jld2_read(SSD_quadratic, rv_gaussian, vb, i)
    rv_symmetric_lsf_array = jld2_read(SSD_quadratic, rv_symmetric_lsf, vb, i)
    rv_asymmetric_lsf_array = jld2_read(SSD_quadratic, rv_asymmetric_lsf, vb, i)

    # rm curve w granulation 
    fig, axs = plt.subplots(2, figsize=(7, 8), sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, line_rv_array, color = 'k', marker = "x", s = 25) 
    rms = round(np.sqrt((np.nansum((line_rv_array - rv_gaussian_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, rv_gaussian_array, color = 'b', linewidth = 2, label = "Gaussian Convolution {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - rv_symmetric_lsf_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, rv_symmetric_lsf_array, color = 'g', linewidth = 2, label = "Symmetric Convolution {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - rv_asymmetric_lsf_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, rv_asymmetric_lsf_array, color = 'r', linewidth = 2, label = "Asymmetric Convolution {}".format(rms))

    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=15)
    axs[0].set_ylabel("RV (m/s)", fontsize=15)
    axs[0].legend(fontsize=11, frameon=False)
    axs[0].tick_params(axis='y', labelsize=11)
    axs[0].set_title("NEID {} $\AA$ RVs".format(line_list[i]))

    # residuals 
    axs[1].scatter(UTC_time, line_rv_array - rv_gaussian_array, color = 'b', s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - rv_symmetric_lsf_array, color = 'g', marker = "x", s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - rv_asymmetric_lsf_array, color = 'r', s = 8) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=15)
    axs[1].set_ylabel("Residuals", fontsize=15) 
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.tick_params(axis='x', which='both', top=True, labeltop=False)
    plt.savefig("rm_{}.png".format(i))
    plt.clf()