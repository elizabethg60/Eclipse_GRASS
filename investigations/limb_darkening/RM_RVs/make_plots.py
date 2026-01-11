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

#these projected models are for N=197, sub=40x40, no extinction (LD only), no sunspots, no MF/CB
SSD_quadratic = h5py.File("projected_SSD_quadratic_gpu.jld2", "r")
RV_SSD_quadratic = SSD_quadratic["RV_list_no_cb"][()]
intensity_SSD_quadratic = SSD_quadratic["intensity"][()]
SSD_4parameter = h5py.File("projected_SSD_4parameter_gpu.jld2", "r")
RV_SSD_4parameter = SSD_4parameter["RV_list_no_cb"][()]
intensity_SSD_4parameter = SSD_4parameter["intensity"][()]

HD_quadratic = h5py.File("projected_HD_gpu_HD_quadratic.jld2", "r")
RV_HD_quadratic = HD_quadratic["RV_list_no_cb"][()]
intensity_HD_quadratic = HD_quadratic["intensity"][()]
HD_4parameter = h5py.File("projected_HD_gpu_HD_4parameter.jld2", "r")
RV_HD_4parameter = HD_4parameter["RV_list_no_cb"][()]
intensity_HD_4parameter = HD_4parameter["intensity"][()]

ld300_quadratic = h5py.File("projected_300_gpu_300_quadratic.jld2", "r")
RV_300_quadratic = ld300_quadratic["RV_list_no_cb"][()]
intensity_300_quadratic = ld300_quadratic["intensity"][()]
ld300_4parameter = h5py.File("projected_300_gpu_300_4parameter.jld2", "r")
RV_300_4parameter = ld300_4parameter["RV_list_no_cb"][()]
intensity_300_4parameter = ld300_4parameter["intensity"][()]

claret_4parameter = h5py.File("projected_SSD_gpu_claret_4parameter.jld2", "r")
RV_claret_4parameter = claret_4parameter["RV_list_no_cb"][()]
intensity_claret_4parameter = claret_4parameter["intensity"][()]

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
    if i == 6:
        continue  
    line_rv_array = jld2_read(line_data, line_rv, vb, i)
    RV_SSD_quadratic_array = jld2_read(SSD_quadratic, RV_SSD_quadratic, vb, i)
    RV_SSD_4parameter_array = jld2_read(SSD_4parameter, RV_SSD_4parameter, vb, i)
    RV_HD_quadratic_array = jld2_read(HD_quadratic, RV_HD_quadratic, vb, i)
    RV_HD_4parameter_array = jld2_read(HD_4parameter, RV_HD_4parameter, vb, i)
    RV_300_quadratic_array = jld2_read(ld300_quadratic, RV_300_quadratic, vb, i)
    RV_300_4parameter_array = jld2_read(ld300_4parameter, RV_300_4parameter, vb, i)
    RV_claret_4parameter_array = jld2_read(claret_4parameter, RV_claret_4parameter, vb, i)

    # rm curve w granulation 
    fig, axs = plt.subplots(2, figsize=(7, 8), sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, line_rv_array, color = 'k', marker = "x", s = 25) 
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_SSD_quadratic_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_SSD_quadratic_array, color = 'b', linewidth = 2, label = "SSD Quadratic {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_SSD_4parameter_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_SSD_4parameter_array, color = 'b', linestyle = "--", linewidth = 2, label = "SSD 4 Parameter {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_HD_quadratic_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_HD_quadratic_array, color = 'r', linewidth = 2, label = "HD Quadratic {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_HD_4parameter_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_HD_4parameter_array, color = 'r', linestyle = "--", linewidth = 2, label = "HD 4 Parameter {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_300_quadratic_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_300_quadratic_array, color = 'g', linewidth = 2, label = "300 Quadratic {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_300_4parameter_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_300_4parameter_array, color = 'g', linestyle = "--", linewidth = 2, label = "300 4 Parameter {}".format(rms))
    rms = round(np.sqrt((np.nansum((line_rv_array - RV_claret_4parameter_array)**2))/len(line_rv_array)),2)
    axs[0].plot(UTC_time, RV_claret_4parameter_array, color = 'purple', linestyle = "--", linewidth = 2, label = "Claret 4 Parameter {}".format(rms))


    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=15)
    axs[0].set_ylabel("RV (m/s)", fontsize=15)
    axs[0].legend(fontsize=11, frameon=False)
    axs[0].tick_params(axis='y', labelsize=11)
    axs[0].set_title("NEID {} $\AA$ RVs".format(line_list[i]))

    # residuals 
    axs[1].scatter(UTC_time, line_rv_array - RV_SSD_quadratic_array, color = 'b', s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_SSD_4parameter_array, color = 'b', marker = "x", s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_HD_quadratic_array, color = 'r', s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_HD_4parameter_array, color = 'r', marker = "x", s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_300_quadratic_array, color = 'g', s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_300_4parameter_array, color = 'g', marker = "x", s = 8) 
    axs[1].scatter(UTC_time, line_rv_array - RV_claret_4parameter_array, color = 'purple', marker = "x", s = 8) 
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=15)
    axs[1].set_ylabel("Residuals", fontsize=15) 
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.tick_params(axis='x', which='both', top=True, labeltop=False)
    plt.savefig("rm_{}.png".format(i))
    plt.clf()