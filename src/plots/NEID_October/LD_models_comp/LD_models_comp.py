import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#GRASS_rv: RV calculated from line profiles + granulation affects on
#GRASS_no_cb: RV calculated from line profiles + granulation affects off

data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

#NL94
#projected rv model - regular
file_regular_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/model_data.jld2", "r")
RV_list_no_cb_NL94 = file_regular_NL94["RV_list_no_cb"][()]
RV_list_cb_NL94  = file_regular_NL94["RV_list_cb"][()]
#projected rv model - quadratic fit
file_quadratic_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/model_data_NL94_quadrant_LD.jld2", "r")
RV_list_no_cb_quadratic_NL94 = file_quadratic_NL94["RV_list_no_cb"][()]
RV_list_cb_quadratic_NL94 = file_quadratic_NL94["RV_list_cb"][()]
#GRASS CB
grass_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_regular.jld2", "r")
lines_NL94 = grass_data_NL94["name"][()]
GRASS_rv_NL94  = grass_data_NL94["rv"][()]
#GRASS no CB
grass_data_no_cb_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_all_lines_rv_off.jld2", "r")
lines_no_cb_NL94 = grass_data_no_cb_NL94["name"][()]
GRASS_no_cb_v_NL94  = grass_data_no_cb_NL94["rv"][()]
#line by line data
line_data_NL94 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/NL94_LD/data/neid_RVlinebyline.jld2", "r")
line_lines_NL94 = line_data_NL94["name"][()]
line_rv_NL94  = line_data_NL94["rv"][()]

#K300
#projected rv model - regular
file_regular_300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/300/data/model_data_Kostogryz_LD_300.jld2", "r")
RV_list_no_cb_300 = file_regular_300["RV_list_no_cb"][()]
RV_list_cb_300 = file_regular_300["RV_list_cb"][()]
#GRASS CB
grass_data_K300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/300/data/neid_all_lines_rv_regular_Kostogryz_LD_300.jld2", "r")
lines_K300 = grass_data_K300["name"][()]
GRASS_rv_K300  = grass_data_K300["rv"][()]
#GRASS no CB
grass_data_no_cb_K300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/300/data/neid_all_lines_rv_off_Kostogryz_LD_300.jld2", "r")
lines_no_cb_K300 = grass_data_no_cb_K300["name"][()]
GRASS_no_cb_v_K300  = grass_data_no_cb_K300["rv"][()]
#line by line data
line_data_K300 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/300/data/neid_RVlinebyline_Kostogryz_LD_300.jld2", "r")
line_lines_K300 = line_data_K300["name"][()]
line_rv_K300  = line_data_K300["rv"][()]

#KSSD
#projected rv model - regular
file_regular_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/SSD/data/model_data_Kostogryz_LD_SSD.jld2", "r")
RV_list_no_cb_SSD = file_regular_SSD["RV_list_no_cb"][()]
RV_list_cb_SSD  = file_regular_SSD["RV_list_cb"][()]
#projected rv model - quadratic fit
file_quadratic_SSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/SSD/data/model_data_Kostogryz_quadrant_LD_SSD.jld2", "r")
RV_list_no_cb_quadratic_SSD = file_quadratic_SSD["RV_list_no_cb"][()]
RV_list_cb_quadratic_SSD = file_quadratic_SSD["RV_list_cb"][()]
#GRASS CB
grass_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/SSD/data/neid_all_lines_rv_regular_Kostogryz_LD_SSD.jld2", "r")
lines_KSSD = grass_data_KSSD["name"][()]
GRASS_rv_KSSD  = grass_data_KSSD["rv"][()]
#GRASS no CB
grass_data_no_cb_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/SSD/data/neid_all_lines_rv_off_Kostogryz_LD_SSD.jld2", "r")
lines_no_cb_KSSD = grass_data_no_cb_KSSD["name"][()]
GRASS_no_cb_v_KSSD  = grass_data_no_cb_KSSD["rv"][()]
#line by line data
line_data_KSSD = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Kostogryz_LD/SSD/data/neid_RVlinebyline_Kostogryz_LD_SSD.jld2", "r")
line_lines_KSSD = line_data_KSSD["name"][()]
line_rv_KSSD  = line_data_KSSD["rv"][()]

RMS_no_CB_NL94 = []
RMS_CB_NL94 = []
RMS_CB_NL94_latter_half = []
RMS_no_CB_K300 = []
RMS_CB_K300 = []
RMS_CB_K300_latter_half = []
RMS_no_CB_KSSD = []
RMS_CB_KSSD = []
RMS_CB_KSSD_latter_half = []

RMS_no_cb_NL94_model = []
RMS_cb_NL94_model = []
RMS_no_cb_quadratic_NL94_model = []
RMS_cb_quadratic_NL94_model = []
RMS_no_cb_300_model = []
RMS_cb_300_model = []
RMS_no_cb_SSD_model = []
RMS_cb_SSD_model  = []
RMS_no_cb_quadratic_SSD_model = []
RMS_cb_quadratic_SSD_model = []
for i in range(0,len(lines_NL94)):
    RV_list_no_cb_NL94_ar = file_regular_NL94[RV_list_no_cb_NL94[i]][()][0:-25]
    RV_list_no_cb_NL94_ar = np.array(RV_list_no_cb_NL94_ar + vb)
    RV_list_no_cb_NL94_ar -= RV_list_no_cb_NL94_ar[-1]

    RV_list_cb_NL94_arr = file_regular_NL94[RV_list_cb_NL94[i]][()][0:-25]
    RV_list_cb_NL94_arr = np.array(RV_list_cb_NL94_arr + vb)
    RV_list_cb_NL94_arr -= RV_list_cb_NL94_arr[-1]

    RV_list_no_cb_quadratic_NL94_ar = file_quadratic_NL94[RV_list_no_cb_quadratic_NL94[i]][()][0:-25]
    RV_list_no_cb_quadratic_NL94_ar = np.array(RV_list_no_cb_quadratic_NL94_ar + vb)
    RV_list_no_cb_quadratic_NL94_ar -= RV_list_no_cb_quadratic_NL94_ar[-1]

    RV_list_cb_quadratic_NL94_ar = file_quadratic_NL94[RV_list_cb_quadratic_NL94[i]][()][0:-25]
    RV_list_cb_quadratic_NL94_ar = np.array(RV_list_cb_quadratic_NL94_ar + vb)
    RV_list_cb_quadratic_NL94_ar -= RV_list_cb_quadratic_NL94_ar[-1]

    RV_list_no_cb_300_ar = file_regular_300[RV_list_no_cb_300[i]][()][0:-25]
    RV_list_no_cb_300_ar = np.array(RV_list_no_cb_300_ar + vb)
    RV_list_no_cb_300_ar -= RV_list_no_cb_300_ar[-1]

    RV_list_cb_300_ar = file_regular_300[RV_list_cb_300[i]][()][0:-25]
    RV_list_cb_300_ar = np.array(RV_list_cb_300_ar + vb)
    RV_list_cb_300_ar -= RV_list_cb_300_ar[-1]

    RV_list_no_cb_SSD_ar = file_regular_SSD[RV_list_no_cb_SSD[i]][()][0:-25]
    RV_list_no_cb_SSD_ar = np.array(RV_list_no_cb_SSD_ar + vb)
    RV_list_no_cb_SSD_ar -= RV_list_no_cb_SSD_ar[-1]

    RV_list_cb_SSD_ar = file_quadratic_NL94[RV_list_cb_SSD[i]][()][0:-25]
    RV_list_cb_SSD_ar = np.array(RV_list_cb_SSD_ar + vb)
    RV_list_cb_SSD_ar -= RV_list_cb_SSD_ar[-1]

    RV_list_no_cb_quadratic_SSD_ar = file_quadratic_SSD[RV_list_no_cb_quadratic_SSD[i]][()][0:-25]
    RV_list_no_cb_quadratic_SSD_ar = np.array(RV_list_no_cb_quadratic_SSD_ar + vb)
    RV_list_no_cb_quadratic_SSD_ar -= RV_list_no_cb_quadratic_SSD_ar[-1]

    RV_list_cb_quadratic_SSD_ar = file_quadratic_SSD[RV_list_cb_quadratic_SSD[i]][()][0:-25]
    RV_list_cb_quadratic_SSD_ar = np.array(RV_list_cb_quadratic_SSD_ar + vb)
    RV_list_cb_quadratic_SSD_ar -= RV_list_cb_quadratic_SSD_ar[-1]

    GRASS_rv_array_NL94 = grass_data_NL94[GRASS_rv_NL94[i]][()][0:-25]
    GRASS_rv_array_NL94 = np.array(GRASS_rv_array_NL94 + vb)
    GRASS_rv_array_NL94 -= GRASS_rv_array_NL94[-1]

    line_rv_array_NL94 = line_data_NL94[line_rv_NL94[i]][()][0:-25]
    line_rv_array_NL94 = np.array(line_rv_array_NL94 + vb)
    line_rv_array_NL94 -= line_rv_array_NL94[-1]

    GRASS_no_cb_array_NL94 = grass_data_no_cb_NL94[GRASS_no_cb_v_NL94[i]][()][0:-25]
    GRASS_no_cb_array_NL94 = np.array(GRASS_no_cb_array_NL94 + vb)
    GRASS_no_cb_array_NL94 -= GRASS_no_cb_array_NL94[-1]

    RMS_no_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_no_cb_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_no_cb_array_NL94))))
    RMS_CB_NL94.append((np.sqrt((np.nansum((line_rv_array_NL94 - GRASS_rv_array_NL94)**2))/len(line_rv_array_NL94 - GRASS_rv_array_NL94))))
    RMS_CB_NL94_latter_half.append((np.sqrt((np.nansum((line_rv_array_NL94[46:-1] - GRASS_rv_array_NL94[46:-1])**2))/len(line_rv_array_NL94[46:-1] - GRASS_rv_array_NL94[46:-1]))))

    RMS_no_cb_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_no_cb_NL94_ar)**2))/len(line_rv_array_NL94 - RV_list_no_cb_NL94_ar))))
    RMS_cb_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_cb_NL94_arr)**2))/len(line_rv_array_NL94 - RV_list_cb_NL94_arr))))
    RMS_no_cb_quadratic_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_no_cb_quadratic_NL94_ar)**2))/len(line_rv_array_NL94 - RV_list_no_cb_quadratic_NL94_ar))))
    RMS_cb_quadratic_NL94_model.append((np.sqrt((np.nansum((line_rv_array_NL94 - RV_list_cb_quadratic_NL94_ar)**2))/len(line_rv_array_NL94 - RV_list_cb_quadratic_NL94_ar))))

    GRASS_rv_array_K300 = grass_data_K300[GRASS_rv_K300[i]][()][0:-25]
    GRASS_rv_array_K300 = np.array(GRASS_rv_array_K300 + vb)
    GRASS_rv_array_K300 -= GRASS_rv_array_K300[-1]

    line_rv_array_K300 = line_data_K300[line_rv_K300[i]][()][0:-25]
    line_rv_array_K300 = np.array(line_rv_array_K300 + vb)
    line_rv_array_K300 -= line_rv_array_K300[-1]

    GRASS_no_cb_array_K300 = grass_data_no_cb_K300[GRASS_no_cb_v_K300[i]][()][0:-25]
    GRASS_no_cb_array_K300 = np.array(GRASS_no_cb_array_K300 + vb)
    GRASS_no_cb_array_K300 -= GRASS_no_cb_array_K300[-1]

    RMS_no_CB_K300.append((np.sqrt((np.nansum((line_rv_array_K300 - GRASS_no_cb_array_K300)**2))/len(line_rv_array_K300 - GRASS_no_cb_array_K300))))
    RMS_CB_K300.append((np.sqrt((np.nansum((line_rv_array_K300 - GRASS_rv_array_K300)**2))/len(line_rv_array_K300 - GRASS_rv_array_K300))))
    RMS_CB_K300_latter_half.append((np.sqrt((np.nansum((line_rv_array_K300[46:-1] - GRASS_rv_array_K300[46:-1])**2))/len(line_rv_array_K300[46:-1] - GRASS_rv_array_K300[46:-1]))))

    RMS_no_cb_300_model.append((np.sqrt((np.nansum((line_rv_array_K300 - RV_list_no_cb_300_ar)**2))/len(line_rv_array_K300 - RV_list_no_cb_300_ar))))
    RMS_cb_300_model.append((np.sqrt((np.nansum((line_rv_array_K300 - RV_list_cb_300_ar)**2))/len(line_rv_array_K300 - RV_list_cb_300_ar))))
    
    GRASS_rv_array_KSSD = grass_data_KSSD[GRASS_rv_KSSD[i]][()][0:-25]
    GRASS_rv_array_KSSD= np.array(GRASS_rv_array_KSSD + vb)
    GRASS_rv_array_KSSD -= GRASS_rv_array_KSSD[-1]

    line_rv_array_KSSD = line_data_KSSD[line_rv_KSSD[i]][()][0:-25]
    line_rv_array_KSSD = np.array(line_rv_array_KSSD + vb)
    line_rv_array_KSSD -= line_rv_array_KSSD[-1]

    GRASS_no_cb_array_KSSD = grass_data_no_cb_KSSD[GRASS_no_cb_v_KSSD[i]][()][0:-25]
    GRASS_no_cb_array_KSSD = np.array(GRASS_no_cb_array_KSSD+ vb)
    GRASS_no_cb_array_KSSD -= GRASS_no_cb_array_KSSD[-1]

    RMS_no_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD - GRASS_no_cb_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_no_cb_array_KSSD))))
    RMS_CB_KSSD.append((np.sqrt((np.nansum((line_rv_array_KSSD- GRASS_rv_array_KSSD)**2))/len(line_rv_array_KSSD - GRASS_rv_array_KSSD))))
    RMS_CB_KSSD_latter_half.append((np.sqrt((np.nansum((line_rv_array_KSSD[46:-1] - GRASS_rv_array_KSSD[46:-1])**2))/len(line_rv_array_KSSD[46:-1] - GRASS_rv_array_KSSD[46:-1]))))
    
    RMS_no_cb_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_no_cb_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_no_cb_SSD_ar))))
    RMS_cb_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_cb_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_cb_SSD_ar))))
    RMS_no_cb_quadratic_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_no_cb_quadratic_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_no_cb_quadratic_SSD_ar))))
    RMS_cb_quadratic_SSD_model.append((np.sqrt((np.nansum((line_rv_array_KSSD - RV_list_cb_quadratic_SSD_ar)**2))/len(line_rv_array_KSSD - RV_list_cb_quadratic_SSD_ar))))


plt.figure(figsize=(12, 6))
plt.scatter(line_lines_NL94, RMS_no_CB_NL94, label = 'CCF RV - no cb')
plt.scatter(line_lines_NL94, RMS_CB_NL94, label = 'CCF RV - GRASS cb')
plt.scatter(line_lines_NL94, RMS_no_cb_NL94_model, label = 'Projected RV - no cb')
plt.scatter(line_lines_NL94, RMS_cb_NL94_model, label = 'Projected RV - Reiners cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("NL94.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_K300, RMS_no_CB_K300, label = 'CCF RV - no cb')
plt.scatter(line_lines_K300, RMS_CB_K300, label = 'CCF RV - GRASS cb')
plt.scatter(line_lines_K300, RMS_no_cb_300_model, label = 'Projected RV - no cb')
plt.scatter(line_lines_K300, RMS_cb_300_model, label = 'Projected RV - Reiners cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("K300.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_KSSD, RMS_no_CB_KSSD, label = 'CCF RV - no cb')
plt.scatter(line_lines_KSSD, RMS_CB_KSSD, label = 'CCF RV - GRASS cb')
plt.scatter(line_lines_KSSD, RMS_no_cb_SSD_model, label = 'Projected RV - no cb')
plt.scatter(line_lines_KSSD, RMS_cb_SSD_model, label = 'Projected RV - Reiners cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("KSSD.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_KSSD, RMS_no_cb_NL94_model, label = 'Projected RV (NL94) - no cb')
plt.scatter(line_lines_KSSD, RMS_cb_NL94_model, label = 'Projected RV (NL94) - GRASS cb')
plt.scatter(line_lines_KSSD, RMS_no_cb_quadratic_NL94_model, label = 'Projected RV (quadratic) - no cb')
plt.scatter(line_lines_KSSD, RMS_cb_quadratic_NL94_model, label = 'Projected RV (quadratic) - Reiners cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("Projected_RV_quadratic_LD_comp_NL94.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_KSSD, RMS_no_cb_SSD_model, label = 'Projected RV (SSD) - no cb')
plt.scatter(line_lines_KSSD, RMS_cb_SSD_model, label = 'Projected RV (SSD) - GRASS cb')
plt.scatter(line_lines_KSSD, RMS_no_cb_quadratic_SSD_model, label = 'Projected RV (quadratic) - no cb')
plt.scatter(line_lines_KSSD, RMS_cb_quadratic_SSD_model, label = 'Projected RV (quadratic) - Reiners cb')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("Projected_RV_quadratic_LD_comp_SSD.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_NL94, RMS_CB_NL94, label = 'NL94')
plt.scatter(line_lines_K300, RMS_CB_K300, label = 'K 300')
plt.scatter(line_lines_KSSD, RMS_CB_KSSD, label = 'K SSD')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("cb_comp.png", bbox_inches='tight')
plt.clf()

plt.figure(figsize=(12, 6))
plt.scatter(line_lines_NL94, RMS_cb_quadratic_NL94_model, label = 'NL94')
plt.scatter(line_lines_KSSD, RMS_cb_quadratic_SSD_model, label = 'SSD')
plt.xlabel("Line Wavelength (Å)", fontsize=12)
plt.ylabel("RV RMS (m/s)", fontsize=12)
plt.xticks(rotation=60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.savefig("cb_comp_quadratic.png", bbox_inches='tight')
plt.clf()