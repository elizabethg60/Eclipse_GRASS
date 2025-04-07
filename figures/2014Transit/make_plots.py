import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime
from astropy.time import Time
from barycorrpy import get_BC_vel

grass_data_combined = h5py.File("projected_SSD_gpu_2014_combined.jld2", "r")
GRASS_rv_combined  = grass_data_combined["RV_list_no_cb"][()]
time_extended  = grass_data_combined["time"][()]
phase_angle  = grass_data_combined["phase_angle"][()]

grass_data_CB = h5py.File("projected_SSD_gpu_2014_CB.jld2", "r")
GRASS_rv_CB  = grass_data_CB["RV_list_no_cb"][()]

grass_data_SH = h5py.File("projected_SSD_gpu_2014_SH.jld2", "r")
GRASS_rv_SH  = grass_data_SH["RV_list_no_cb"][()]

model1 = h5py.File("projected_SSD_gpu_2014_model1_275.jld2", "r")
model1_rv  = model1["RV_list_no_cb"][()]

UTC_time_ex = []
time_julian_ext = []
for i in time_extended:
    dt = datetime.strptime(i.decode("utf-8"), "%Y-%m-%dT%H:%M:%S.%f")
    UTC_time_ex.append(dt)
    time_julian_ext.append((Time(dt)).jd)

vb_ext, warnings, flag = get_BC_vel(JDUTC=time_julian_ext, lat=28.754081 , longi=-17.889242, alt=235.9, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

def jld2_read(jld2_file, variable, v, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + v)
    array -= array[-1]
    return array

df = pd.read_csv('data.csv')
t = Time(df['x'], format='mjd')
# Convert to a regular datetime object
datetime_obj = t.to_datetime()
data = df['y']

for i in range(0, 1):
    
    model1_arr = jld2_read(model1, model1_rv, vb_ext, i)
    GRASS_rv_array_SH = jld2_read(grass_data_SH, GRASS_rv_SH, vb_ext, i)
    GRASS_rv_array_CB = jld2_read(grass_data_CB, GRASS_rv_CB, vb_ext, i)
    GRASS_rv_array_combined = jld2_read(grass_data_combined, GRASS_rv_combined, vb_ext, i)
    phase_angle_array = grass_data_combined[phase_angle[i]][()]
    differences = np.abs(phase_angle_array - 0.4)
    indices = np.argsort(differences)[:2]

    # rm curve 
    fig, ax = plt.subplots()
    ax.plot(UTC_time_ex[65:-45], GRASS_rv_array_SH[65:-45], color = 'b', linewidth = 3, label = "SH")
    # ax.plot(UTC_time_ex[65:-45], GRASS_rv_array_combined[65:-45], color = 'r', linewidth = 3, label = "Combined")
    ax.plot(UTC_time_ex[65:-45], GRASS_rv_array_CB[65:-45], color = 'r', linewidth = 3, label = "CB")
    # ax.plot(UTC_time_ex[65:-45], model1_arr[65:-45], color = 'k', linewidth = 1, label = "Model 1")
    ax.scatter(datetime_obj, data, color = 'gray')
    # start_date = pd.Timestamp('2014-01-05T15:38:03')
    # end_date = pd.Timestamp('2014-01-06T02:00:00')
    # ax.axvspan(start_date, end_date, color='gray', alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel("Time on 01/05/2014 - 01/07/2014 (UTC)") #full 01/04/2014 - 01/08/2014
    ax.set_ylabel("RV [m/s]")
    ax.legend()
    plt.savefig("phase_curve_models_proposal.png")
    plt.clf()

    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].plot(UTC_time_ex[65:-45], GRASS_rv_array_combined[65:-45], color = 'r', linewidth = 2, label = "Combined")
    # axs[0].plot(UTC_time_ex[65:-45], GRASS_rv_array_SH[65:-45], color = 'b', linewidth = 2, label = "SH")
    # axs[0].plot(UTC_time_ex[65:-45], GRASS_rv_array_CB[65:-45], color = 'g', linewidth = 2, label = "CB")
    # axs[0].scatter(datetime_obj, data, color = 'k')
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # start_date = pd.Timestamp('2014-01-05T15:38:03')
    # end_date = pd.Timestamp('2014-01-06T02:00:00')
    # axs[0].axvspan(start_date, end_date, color='gray', alpha=0.3)
    # axs[1].axvspan(start_date, end_date, color='gray', alpha=0.3)
    # axs[0].set_xlabel("Time on 01/05/2014 - 01/07/2014 (UTC)", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time_ex[65:-45], phase_angle_array[65:-45], color = 'r', marker = "x", s = 3)  
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("Time on 01/05/2014 - 01/07/2014 (UTC)", fontsize=12)
    # axs[1].set_ylabel("Phase Angle (degrees)", fontsize=12) 
    # # axs[0].axvline(x = UTC_time_ex[indices[0]])
    # # axs[0].axvline(x = UTC_time_ex[indices[1]])
    # # axs[1].axvline(x = UTC_time_ex[indices[0]])
    # # axs[1].axvline(x = UTC_time_ex[indices[1]])
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("phase_angles_zoom")
    # plt.clf()