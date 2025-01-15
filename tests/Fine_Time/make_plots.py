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
grass_data_no_cb = h5py.File("neid_all_lines_rv_off_SSD_gpu_fine_time.jld2", "r")
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

def jld2_read_fine_time(jld2_file, variable, vb, index):
    original = jld2_file[variable[index]][()]
    original = np.array(original + vb_fine)
    original -= original[-1]

    array = [(original[i:i+56]) for i in range(0, len(original), 56)]
    for i in range(0, len(array)):
        plt.scatter(range(0, 56), array[i])
        plt.savefig("1sec_RV/RV_single_exposure_{}".format(i))
        plt.clf()

def jld2_read_brightness(jld2_file, variable, index):
    original = jld2_file[variable[index]][()]
    array = [(original[i:i+56]) for i in range(0, len(original), 56)]
    for i in range(0, len(array)):
        plt.scatter(range(0, 56), array[i])
        plt.savefig("1sec_brightness/brightness_{}".format(i))
        plt.clf()

jld2_read_fine_time(grass_data_no_cb, GRASS_no_cb, vb_fine, 0)
brightness_array = jld2_read_brightness(grass_data_no_cb, brightness, 0)
