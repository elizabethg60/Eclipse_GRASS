import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel
from scipy.interpolate import CubicSpline
import random

grass_data_combined = h5py.File("projected_SSD_gpu_2026_combined.jld2", "r")
GRASS_rv_combined  = grass_data_combined["RV_list_no_cb"][()]
time_extended  = grass_data_combined["time"][()]
phase_angle  = grass_data_combined["phase_angle"][()]

grass_data_CB = h5py.File("projected_SSD_gpu_2026_CB.jld2", "r")
GRASS_rv_CB  = grass_data_CB["RV_list_no_cb"][()]

grass_data_SH = h5py.File("projected_SSD_gpu_2026_SH.jld2", "r")
GRASS_rv_SH  = grass_data_SH["RV_list_no_cb"][()]

UTC_time_ex = []
time_julian_ext = []
for i in time_extended:
    dt = datetime.strptime(i.decode("utf-8"), "%Y-%m-%dT%H:%M:%S.%f")
    UTC_time_ex.append(dt)
    time_julian_ext.append((Time(dt)).jd)

vb_ext, warnings, flag = get_BC_vel(JDUTC=time_julian_ext, lat=31.9583, longi= -111.5967, alt=2.097938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

def jld2_read(jld2_file, variable, v, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + v)
    array -= array[-1]
    return array

df = pd.read_csv(str("2026Transit.txt"), header=None)

for i in range(0, 1):
    
    GRASS_rv_array_SH = jld2_read(grass_data_SH, GRASS_rv_SH, vb_ext, i)
    GRASS_rv_array_CB = jld2_read(grass_data_CB, GRASS_rv_CB, vb_ext, i)
    GRASS_rv_array_combined = jld2_read(grass_data_combined, GRASS_rv_combined, vb_ext, i)
    phase_angle_array = grass_data_combined[phase_angle[i]][()]
    differences = np.abs(phase_angle_array - 0.4)
    indices = np.argsort(differences)[:2]

    # # rm curve 
    # fig, ax = plt.subplots()
    # ax.plot(UTC_time_ex[67:-47], GRASS_rv_array_SH[67:-47], color = 'b', linewidth = 3, label = "SH")
    # # ax.plot(UTC_time_ex[67:-47], GRASS_rv_array_combined[67:-47], color = 'r', linewidth = 3, label = "Combined")
    # ax.plot(UTC_time_ex[67:-47], GRASS_rv_array_CB[67:-47], color = 'r', linewidth = 3, label = "CB")
    # # start_date = pd.Timestamp('2026-01-10 07:00:23')
    # # end_date = pd.Timestamp('2026-01-10 11:45:56')
    # # ax.axvspan(start_date, end_date, color='gray', alpha=0.3) #transit
    # #visibile from NEID
    # #gray - night time in AZ (2-13 UTC) & airmass < 2.5
    # start_date = pd.Timestamp('2026-01-09 06:00')
    # end_date = pd.Timestamp('2026-01-09 12:30')
    # ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    # start_date = pd.Timestamp('2026-01-10 02:30')
    # end_date = pd.Timestamp('2026-01-10 12:30')
    # ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    # start_date = pd.Timestamp('2026-01-11 02:30')
    # end_date = pd.Timestamp('2026-01-11 12:30')
    # ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    # # #green - Jupiter not an issue
    # # start_date = pd.Timestamp('2026-01-09 06:00')
    # # end_date = pd.Timestamp('2026-01-09 18:00')
    # # ax.axvspan(start_date, end_date, color='green', alpha=0.4)
    # # start_date = pd.Timestamp('2026-01-10 08:00')
    # # end_date = pd.Timestamp('2026-01-11 13:00')
    # # ax.axvspan(start_date, end_date, color='green', alpha=0.4)
    # ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # ax.set_xlabel("Time on 01/09/2026 - 01/11/2026 (UTC)") 
    # ax.set_ylabel("RV [m/s]")
    # ax.legend()
    # plt.savefig("phase_curve_models_proposal_vis.png")
    # plt.clf()

# Define start and end datetimes
    dates = UTC_time_ex[67:-47]
    start_datetime = datetime(2026, 1, 9, 6, 0, 0)
    end_datetime = datetime(2026, 1, 9, 12, 30, 0)
    time_step = timedelta(minutes=5.4)
    # Generate an array of datetime objects from start to end with the given time step
    current_datetime = start_datetime
    datetimes1 = []
    datestimes_datetime1 = []

    while current_datetime <= end_datetime:
        datestimes_datetime1.append(current_datetime + time_step)
        datetimes1.append(((current_datetime + time_step) - dates[0]).total_seconds())
        current_datetime += time_step

    start_datetime = datetime(2026, 1, 10, 2, 30, 0)
    end_datetime = datetime(2026, 1, 10, 12, 30, 0)
    current_datetime = start_datetime
    datetimes2 = []
    datestimes_datetime2 = []

    while current_datetime <= end_datetime:
        datestimes_datetime2.append(current_datetime + time_step)
        datetimes2.append(((current_datetime + time_step) - dates[0]).total_seconds())
        current_datetime += time_step

    start_datetime = datetime(2026, 1, 11, 2, 30, 0)
    end_datetime = datetime(2026, 1, 11, 12, 30, 0)
    current_datetime = start_datetime
    datetimes3 = []
    datestimes_datetime3 = []

    while current_datetime <= end_datetime:
        datestimes_datetime3.append(current_datetime + time_step)
        datetimes3.append(((current_datetime + time_step) - dates[0]).total_seconds())
        current_datetime += time_step

    # Convert datetime to numerical values (e.g., number of days since the first date)
    date_nums = []
    for i in dates:
        date_nums.append((i - dates[0]).total_seconds())

    # Create the cubic spline interpolator
    cs = CubicSpline(date_nums, GRASS_rv_array_SH[67:-47] - GRASS_rv_array_CB[67:-47])

    def estimate_rv(timestamps, timevalue):
        x_values = cs(timevalue)
        new_time = timestamps[::10]
        random_value = []
        for i in range(0, len(new_time)):
            random_value.append(random.gauss(np.mean(x_values), 0.421))
        return new_time, random_value

    

    # rm curve 
    fig, ax = plt.subplots()
    ax.errorbar(estimate_rv(datestimes_datetime1, datetimes1)[0], estimate_rv(datestimes_datetime1, datetimes1)[1], yerr=[0.421]*len(estimate_rv(datestimes_datetime1, datetimes1)[0]), fmt='o', color = 'red')
    ax.errorbar(estimate_rv(datestimes_datetime2, datetimes2)[0], estimate_rv(datestimes_datetime2, datetimes2)[1], yerr=[0.421]*len(estimate_rv(datestimes_datetime2, datetimes2)[0]), fmt='o', color = 'red')
    ax.errorbar(estimate_rv(datestimes_datetime3, datetimes3)[0], estimate_rv(datestimes_datetime3, datetimes3)[1], yerr=[0.421]*len(estimate_rv(datestimes_datetime3, datetimes3)[0]), fmt='o', color = 'red')
    ax.plot(UTC_time_ex[67:-47], GRASS_rv_array_SH[67:-47] - GRASS_rv_array_CB[67:-47], color = 'k', linewidth = 3, label = "SH-CB")
    #visibile from NEID
    #gray - night time in AZ (2-13 UTC) & airmass < 2.5
    start_date = pd.Timestamp('2026-01-09 06:00')
    end_date = pd.Timestamp('2026-01-09 12:30')
    ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    start_date = pd.Timestamp('2026-01-10 02:30')
    end_date = pd.Timestamp('2026-01-10 12:30')
    ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    start_date = pd.Timestamp('2026-01-11 02:30')
    end_date = pd.Timestamp('2026-01-11 12:30')
    ax.axvspan(start_date, end_date, color='gray', alpha=0.4)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel("Time on 01/09/2026 - 01/11/2026 (UTC)") 
    ax.set_ylabel("RV [m/s]")
    ax.legend()
    plt.savefig("phase_curve_models_proposal_residuals.png")
    plt.clf()

    # fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    # axs[0].plot(UTC_time_ex, GRASS_rv_array_combined, color = 'r', linewidth = 2, label = "Combined")
    # axs[0].plot(UTC_time_ex, GRASS_rv_array_SH, color = 'b', linewidth = 2, label = "SH")
    # axs[0].plot(UTC_time_ex, GRASS_rv_array_CB, color = 'g', linewidth = 2, label = "CB")
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # start_date = pd.Timestamp('2026-01-10 07:00:23')
    # end_date = pd.Timestamp('2026-01-10 11:45:56')
    # axs[0].axvspan(start_date, end_date, color='gray', alpha=0.3)
    # axs[1].axvspan(start_date, end_date, color='gray', alpha=0.3)
    # axs[0].set_xlabel("Time on 01/08/2026 - 01/12/2026 (UTC)", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # #residuals
    # axs[1].scatter(UTC_time_ex, phase_angle_array, color = 'r', marker = "x", s = 3)  
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("Time on 01/08/2026 - 01/12/2026 (UTC)", fontsize=12)
    # axs[1].set_ylabel("Phase Angle (degrees)", fontsize=12) 
    # axs[0].axvline(x = UTC_time_ex[indices[0]])
    # axs[0].axvline(x = UTC_time_ex[indices[1]])
    # axs[1].axvline(x = UTC_time_ex[indices[0]])
    # axs[1].axvline(x = UTC_time_ex[indices[1]])
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig("phase_angles")
    # plt.clf()   