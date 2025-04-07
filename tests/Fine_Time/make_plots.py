import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel

# read in data
# GRASS
# grass_data_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_time/SSD_3ext/data/neid_all_lines_rv_off_SSD_gpu_fine_time_first_line_N300_sub40.jld2", "r")
# #grid: N = 300, sub = 40 / for first line 5250 / SSD / every 1 sec / GPU / extinction 
grass_data_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_time_5434_line_N900_sub40.jld2", "r")
#grid: N = 900, sub = 40 / for 5434 / SSD / every 1 sec / GPU / no extinction 
GRASS_no_cb  = grass_data_no_cb["rv"][()]
GRASS_no_cb_var  = grass_data_no_cb["rv_var"][()]
GRASS_no_cb_lsf  = grass_data_no_cb["rv_lsf"][()]
GRASS_no_cb_var_lsf  = grass_data_no_cb["rv_var_lsf"][()]
brightness = grass_data_no_cb["brightness"][()]

# data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")
fine_sample_UTC_time = []
for i in data["obsdate"][15:-178]:
    inner = []
    for sec in np.arange(0.0, 56.0, 1.0):
        dt = (datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=sec)).strftime("%Y-%m-%dT%H:%M:%S")
        inner.append((Time(dt)).jd)
    fine_sample_UTC_time.append(inner)
fine_time_julian = sum(fine_sample_UTC_time, [])
vb_fine, warnings, flag = get_BC_vel(JDUTC=fine_time_julian, lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

def jld2_read_rv(jld2_file, variable, index):
    original = jld2_file[variable[index]][()]
    original = np.array(original + vb_fine)
    original -= original[-1]
    array = [(original[i:i+56]) for i in range(0, len(original), 56)]
    slope = []

    for i in range(0, len(array)):
        plt.scatter(range(0, 56), array[i], color = "g")
        plt.xlabel("Time within 55 s exposure (s)")
        plt.ylabel("RV (m/s)")
        plt.savefig("FeI5434_N900_sub40/NoExt/RV/RV_{}".format(i))
        plt.clf()
        slope.append((array[i][55] - array[i][0])/len(array[i])) 

    plt.scatter(range(0, len(slope)), slope, color = "r")
    plt.xlabel("Count of 55 s exposure")
    plt.ylabel("RV Slope")
    plt.savefig("FeI5434_N900_sub40/NoExt/RV/slope", bbox_inches='tight')
    plt.clf()  

def jld2_read_brightness(jld2_file, variable, index):
    original = jld2_file[variable[index]][()]
    original = original / np.max(original)
    array = [(original[i:i+56]) for i in range(0, len(original), 56)]
    slope = []

    for i in range(0, len(array)):
        plt.scatter(range(0, 56), array[i], color = "g")
        plt.xlabel("Time within 55 s exposure (s)")
        plt.ylabel("Relative Flux")
        plt.savefig("FeI5434_N900_sub40/NoExt/Flux/Flux_{}".format(i))
        plt.clf()
        slope.append((array[i][55] - array[i][0])/len(array[i]))

    plt.scatter(range(0, len(slope)), slope, color = "r")
    plt.xlabel("Count of 55 s exposure")
    plt.ylabel("Flux Slope")
    plt.savefig("FeI5434_N900_sub40/NoExt/Flux/slope", bbox_inches='tight')
    plt.clf()  

def jld2_read_rv_convergence(jld2_file, variable1, variable2, N_range_Ngrid):
    variable1_rms = []
    variable2_rms = []
    for i in range(0,len(N_range_Ngrid)): 
        rv = jld2_file[variable1[i]][()]
        rv = np.array(rv + vb_fine[0:len(rv)])
        slope_rv, intercept_rv = np.polyfit(range(0, len(rv)), rv, 1)
        fit_rv = slope_rv * range(0, len(rv)) + intercept_rv
        plt.scatter(range(0, len(rv)), rv)
        plt.plot(range(0, len(rv)), fit_rv)
        plt.xlabel("Time within 55 s exposure (s)")
        plt.ylabel("RV (m/s)")
        plt.savefig("FeI5434_Nconvergence_sub40_no_ext/RV/N_{}".format(N_range_Ngrid[i]))
        plt.clf()
        variable1_rms.append(np.sqrt((np.nansum((fit_rv - rv)**2))/len(fit_rv)))

        flux = jld2_file[variable2[i]][()]
        slope_flux, intercept_flux = np.polyfit(range(0, len(flux)), flux, 1)
        fit_flux = slope_flux * range(0, len(flux)) + intercept_flux
        plt.scatter(range(0, len(flux)), flux)
        plt.plot(range(0, len(flux)), fit_flux)
        plt.xlabel("Time within 55 s exposure (s)")
        plt.ylabel("Flux")
        plt.savefig("FeI5434_Nconvergence_sub40_no_ext/Flux/N_{}".format(N_range_Ngrid[i]))
        plt.clf()
        variable2_rms.append(np.sqrt((np.nansum((fit_flux - flux)**2))/len(fit_flux)))

    plt.scatter(N_range_Ngrid, variable1_rms, color = "r")
    plt.xlabel("N (latitude slices)")
    plt.ylabel("Residual RV rms (m/s)")
    plt.savefig("FeI5434_Nconvergence_sub40_no_ext/RV/convergence", bbox_inches='tight')
    plt.clf() 

    plt.scatter(N_range_Ngrid, variable2_rms, color = "r")
    plt.xlabel("N (latitude slices)")
    plt.ylabel("Residual Flux rms")
    plt.savefig("FeI5434_Nconvergence_sub40_no_ext/Flux/convergence", bbox_inches='tight')
    plt.clf() 

jld2_read_rv(grass_data_no_cb, GRASS_no_cb, 8)
jld2_read_brightness(grass_data_no_cb, brightness, 0)

# N_convergence = h5py.File("FeI5434_Nconvergence_sub40_no_ext.jld2", "r")
# N_range_Ngrid = np.arange(300,1010,10)
# rv_Ngrid  = N_convergence["rv"][()]
# brightness_Ngrid  = N_convergence["brightness"][()]

# jld2_read_rv_convergence(N_convergence, rv_Ngrid, brightness_Ngrid, N_range_Ngrid)

# note: over lap of flux and RV are the same 