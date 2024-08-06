import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

#NL94
#projected rv model - regular
file_regular_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_october_N_50_NL94.jld2", "r")
RV_list_no_cb_NL94 = file_regular_NL94["RV_list_no_cb"][()]
#GRASS CB
grass_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_regular_NL94.jld2", "r")
lines_NL94 = grass_data_NL94["name"][()]
GRASS_rv_NL94  = grass_data_NL94["rv"][()]
#GRASS no CB
grass_data_no_cb_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_off_NL94.jld2", "r")
lines_no_cb_NL94 = grass_data_no_cb_NL94["name"][()]
GRASS_no_cb_v_NL94  = grass_data_no_cb_NL94["rv"][()]
#line by line data
line_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_RVlinebyline_NL94.jld2", "r")
line_lines_NL94 = line_data_NL94["name"][()]
line_rv_NL94  = line_data_NL94["rv"][()]

#KSSD
#projected rv model - regular
file_regular_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/data/neid_october_N_50_KSSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
#GRASS CB
grass_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/data/neid_all_lines_rv_regular_KSSD.jld2", "r")
lines_KSSD = grass_data_KSSD["name"][()]
GRASS_rv_KSSD  = grass_data_KSSD["rv"][()]
#GRASS no CB
grass_data_no_cb_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/data/neid_all_lines_rv_off_KSSD.jld2", "r")
lines_no_cb_KSSD = grass_data_no_cb_KSSD["name"][()]
GRASS_no_cb_v_KSSD  = grass_data_no_cb_KSSD["rv"][()]
#line by line data
line_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/data/neid_RVlinebyline_Kostogryz_LD_SSD.jld2", "r")
line_lines_KSSD = line_data_KSSD["name"][()]
line_rv_KSSD  = line_data_KSSD["rv"][()]

RMS_no_CB_NL94 = []
RMS_CB_NL94 = []
RMS_no_CB_KSSD = []
RMS_CB_KSSD = []
RMS_no_cb_NL94_model = []
RMS_no_cb_SSD_model = []

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    array -= array[-1]
    return array

for i in range(0,len(lines_NL94)):
    RV_list_no_cb_NL94_ar = jld2_read(file_regular_NL94, RV_list_no_cb_NL94, vb, i)
    RV_list_no_cb_SSD_ar = jld2_read(file_regular_SSD, RV_list_no_cb_SSD, vb, i)
    GRASS_rv_array_NL94 = jld2_read(grass_data_NL94, GRASS_rv_NL94, vb, i)
    line_rv_array_NL94 = jld2_read(line_data_NL94, line_rv_NL94, vb, i)
    GRASS_no_cb_array_NL94 = jld2_read(grass_data_no_cb_NL94, GRASS_no_cb_v_NL94, vb, i)
    GRASS_rv_array_KSSD = jld2_read(grass_data_KSSD, GRASS_rv_KSSD, vb, i)
    line_rv_array_KSSD = jld2_read(line_data_KSSD, line_rv_KSSD, vb, i)
    GRASS_no_cb_array_KSSD = jld2_read(grass_data_no_cb_KSSD, GRASS_no_cb_v_KSSD, vb, i)

    RMS_no_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_no_cb_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_no_cb_array_NL94))))
    RMS_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_rv_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_rv_array_NL94))))
    RMS_no_cb_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_no_cb_NL94_ar)**2))/len(line_rv_array_NL94 - RV_list_no_cb_NL94_ar))))
    RMS_no_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD - GRASS_no_cb_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_no_cb_array_KSSD))))
    RMS_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD- GRASS_rv_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_rv_array_KSSD))))
    RMS_no_cb_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_no_cb_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_no_cb_SSD_ar))))


plt.figure(figsize=(12, 6))
plt.scatter(line_lines_NL94, RMS_no_CB_NL94, label = 'CCF RV - no cb')
plt.scatter(line_lines_NL94, RMS_CB_NL94, label = 'CCF RV - GRASS cb')
plt.scatter(line_lines_NL94, RMS_no_cb_NL94_model, label = 'Projected RV - no cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("NL94.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_KSSD, RMS_no_CB_KSSD, label = 'CCF RV - no cb')
plt.scatter(line_lines_KSSD, RMS_CB_KSSD, label = 'CCF RV - GRASS cb')
plt.scatter(line_lines_KSSD, RMS_no_cb_SSD_model, label = 'Projected RV - no cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("KSSD.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_NL94, RMS_CB_NL94, label = 'NL94')
plt.scatter(line_lines_KSSD, RMS_CB_KSSD, label = 'KSSD')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("cb_comp.png", bbox_inches='tight')
plt.clf()