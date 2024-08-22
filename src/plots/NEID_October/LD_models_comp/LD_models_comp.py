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

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)[0:-25]

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

#HD
#projected rv model - regular
file_regular_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_october_N_50_HD.jld2", "r")
RV_list_no_cb_HD = file_regular_HD["RV_list_no_cb"][()]
#GRASS CB
grass_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_all_lines_rv_regular_HD.jld2", "r")
lines_HD = grass_data_HD["name"][()]
GRASS_rv_HD  = grass_data_HD["rv"][()]
#GRASS no CB
grass_data_no_cb_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_all_lines_rv_off_HD.jld2", "r")
lines_no_cb_HD = grass_data_no_cb_HD["name"][()]
GRASS_no_cb_v_HD  = grass_data_no_cb_HD["rv"][()]
#line by line data
line_data_HD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KHD/data/neid_RVlinebyline_HD.jld2", "r")
line_lines_HD = line_data_HD["name"][()]
line_rv_HD  = line_data_HD["rv"][()]

RMS_no_CB_NL94 = []
RMS_CB_NL94 = []
RMS_no_CB_KSSD = []
RMS_CB_KSSD = []
RMS_no_CB_HD = []
RMS_CB_HD = []
RMS_no_cb_NL94_model = []
RMS_no_cb_SSD_model = []
RMS_no_cb_HD_model = []
RMS_CB_NL94_pipeline = []
RMS_CB_KSSD_pipeline = []
RMS_CB_HD_pipeline = []

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    array -= array[-1]
    return array

lines = []
for i in range(0,len(lines_NL94)):
    if i == 5:
        continue
    lines.append(line_lines_NL94[i])
    RV_list_no_cb_NL94_ar = jld2_read(file_regular_NL94, RV_list_no_cb_NL94, vb, i)
    RV_list_no_cb_SSD_ar = jld2_read(file_regular_SSD, RV_list_no_cb_SSD, vb, i)
    GRASS_rv_array_NL94 = jld2_read(grass_data_NL94, GRASS_rv_NL94, vb, i)
    RV_list_no_cb_HD_ar = jld2_read(file_regular_HD, RV_list_no_cb_HD, vb, i)
    line_rv_array_NL94 = jld2_read(line_data_NL94, line_rv_NL94, vb, i)
    GRASS_no_cb_array_NL94 = jld2_read(grass_data_no_cb_NL94, GRASS_no_cb_v_NL94, vb, i)
    GRASS_rv_array_KSSD = jld2_read(grass_data_KSSD, GRASS_rv_KSSD, vb, i)
    line_rv_array_KSSD = jld2_read(line_data_KSSD, line_rv_KSSD, vb, i)
    GRASS_no_cb_array_KSSD = jld2_read(grass_data_no_cb_KSSD, GRASS_no_cb_v_KSSD, vb, i)
    GRASS_rv_array_HD = jld2_read(grass_data_HD, GRASS_rv_HD, vb, i)
    line_rv_array_HD = jld2_read(line_data_HD, line_rv_HD, vb, i)
    GRASS_no_cb_array_HD = jld2_read(grass_data_no_cb_HD, GRASS_no_cb_v_HD, vb, i)

    RMS_no_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_no_cb_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_no_cb_array_NL94))))
    RMS_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_rv_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_rv_array_NL94))))
    RMS_no_cb_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_no_cb_NL94_ar)**2))/len(line_rv_array_NL94 - RV_list_no_cb_NL94_ar))))
    RMS_no_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD - GRASS_no_cb_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_no_cb_array_KSSD))))
    RMS_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD- GRASS_rv_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_rv_array_KSSD))))
    RMS_no_cb_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_no_cb_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_no_cb_SSD_ar))))
    RMS_no_CB_HD.append((np.sqrt((np.nansum((line_rv_array_HD - GRASS_no_cb_array_HD)**2))/len(line_rv_array_HD - GRASS_no_cb_array_HD))))
    RMS_CB_HD.append((np.sqrt((np.nansum((line_rv_array_HD- GRASS_rv_array_HD)**2))/len(line_rv_array_HD - GRASS_rv_array_HD))))
    RMS_no_cb_HD_model.append((np.sqrt((np.nansum((line_rv_array_HD - RV_list_no_cb_HD_ar)**2))/len(line_rv_array_HD - RV_list_no_cb_HD_ar))))
    RMS_CB_NL94_pipeline.append((np.sqrt((np.nansum((rv_obs - GRASS_rv_array_NL94)**2))/len(rv_obs - GRASS_rv_array_NL94))))
    RMS_CB_KSSD_pipeline.append((np.sqrt((np.nansum((rv_obs - GRASS_rv_array_KSSD)**2))/len(rv_obs - GRASS_rv_array_KSSD))))
    RMS_CB_HD_pipeline.append((np.sqrt((np.nansum((rv_obs - GRASS_rv_array_HD)**2))/len(rv_obs - GRASS_rv_array_HD))))

plt.figure(figsize=(12, 6))
plt.scatter(lines, RMS_no_CB_NL94, label = 'CCF RV - no cb')
plt.scatter(lines, RMS_CB_NL94, label = 'CCF RV - GRASS cb')
plt.scatter(lines, RMS_no_cb_NL94_model, label = 'Projected RV - no cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("NL94.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(lines, RMS_no_CB_KSSD, label = 'CCF RV - no cb')
plt.scatter(lines, RMS_CB_KSSD, label = 'CCF RV - GRASS cb')
plt.scatter(lines, RMS_no_cb_SSD_model, label = 'Projected RV - no cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("KSSD.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(lines, RMS_no_CB_HD, label = 'CCF RV - no cb')
plt.scatter(lines, RMS_CB_HD, label = 'CCF RV - GRASS cb')
plt.scatter(lines, RMS_no_cb_HD_model, label = 'Projected RV - no cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("HD.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(lines, RMS_CB_NL94, label = 'NL94 - line')
plt.scatter(lines, RMS_CB_KSSD, label = 'KSSD - line')
plt.scatter(lines, RMS_CB_HD, label = 'HD - line')
plt.scatter(lines, RMS_CB_NL94_pipeline, label = 'NL94 - pipeline')
plt.scatter(lines, RMS_CB_KSSD_pipeline, label = 'KSSD - pipeline')
plt.scatter(lines, RMS_CB_HD_pipeline, label = 'HD - pipeline')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("cb_comp.png", bbox_inches='tight')
plt.clf()