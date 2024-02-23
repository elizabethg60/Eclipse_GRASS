import h5py
import math
import csv
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt 
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from astropy.time import Time

GRASS_rv = [108.64094225951013, 107.88620410583766, 107.6513559241415, 108.67327613716924, 109.33710070113621, 112.07159149471435, 111.90315663898716, 117.2990576939628, 118.09502309249397, 120.47932540282963, 117.22694558126557, 120.86047032618471, 123.45851372345118, 121.73625038675146, 121.87777603722893, 121.17191741772321, 123.10664889154755, 126.63463844828028, 127.91857842780672, 132.002738195899, 133.60697990106024, 131.11943328112534, 131.5098414633938, 138.09201302996266, 137.81543550091402, 137.83232339062457, 137.52737553653466, 138.94903917058207, 141.81370241913172, 144.44691699606904, 146.01035239693155, 148.18358629884443, 142.8342707551814, 147.51536134722707, 148.07756637720416, 151.42841604573545, 151.92162230576545, 152.70320076523853, 155.5246061144359, 157.18567441493414, 160.1420695133973, 161.09379386101384, 162.9784926771806, 165.2793855810665, 166.3400290779997, 167.0711035758846, 171.03436370102597, 167.7019890561081, 169.1093936760339, 161.89553422834356, 148.69603091963623, 135.97695652525502, 119.96684325564651, 100.31660562749487, 78.82766164917521, 53.312961825234105, 27.67201621955743, 4.652183204235477, -23.44454854835339, -52.51418046073086, -78.0947647174172, -109.09178939762806, -142.3180905761803, -171.4484604062677, -206.77594034304235, -225.80572041527032, -263.9605245580223, -287.5815330573805, -319.61529976460037, -342.45443815491825, -373.5680026083804, -390.6782468709241, -414.0750397794598, -415.9039292328149, -424.4517561879557, -423.95675523449313, -414.20889649454523, -401.3649299046067, -359.15732924425834, -303.24247663383795, -230.49037434504496, -145.7591907766079, -31.050079197115608, 106.24064846284811, 238.6546397792086, 374.1530677504269, 515.1261816348315, 638.2137318470158, 728.6023238115529, 803.1368917649849, 864.5522077169096, 901.2501902989211, 928.0313925183257, 944.474796032202, 953.8264311609048, 942.0500449463013, 935.526508688026, 911.7871947477822, 895.0746113009695, 872.6700675884048, 838.9501254600683, 812.2021445889585, 786.8114722705164, 755.0321109542416, 725.3963419576933, 699.33410457612, 663.3754542810955, 637.5574351864122, 600.6980173915163, 575.1794036732814, 546.2483751690751, 518.7795700502664, 491.09308717793806, 468.5986391685009, 442.7894589820534, 424.3652463995164, 404.322722913144, 383.4006448782058, 368.3896542099512, 352.76594568734157, 343.7299486001338, 330.1415218853817, 324.9752201381821, 325.5058390060303, 324.97470655564297, 329.94115613267417, 336.55494401368014, 335.0970879015002, 334.1618937386518, 341.99027392357084, 340.26078730670366, 343.60913317717325, 341.5571860964697, 348.35770536820456, 352.4882133381904, 354.7788413163064, 352.9623932171206, 356.3273418526258, 361.75778840352297, 363.098660746068, 364.362829756919, 365.62082619060743, 365.31543951483155, 373.82341201716264, 372.9558031833476, 372.65146579754986, 377.8442759872173, 379.4996545531012, 381.8655116388073, 383.7441018412137, 397.69030042890205, 399.7022553136386, 401.50167457810363, 403.9673445616713, 407.8162894404973, 405.6231365240052, 408.85852815996935, 411.5100144352195, 414.70821392937523]

#read in data
#model
#file = h5py.File("Resolution/model_data100_20sub.jld2", "r")
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
RV_list_cb_new  = file["RV_list_cb_new"][()]
intensity_list = file["intensity_list"][()]
vel_no_cb = file["vel_no_cb"][()]
vel_cb = file["vel_cb"][()]

#data 
f=open("Reiners_Data.txt","r")
lines=f.readlines()[1:]
time=[]
raw_rv = []
time_julian = []
UTC_time = []
for x in lines:
    dt = datetime.strptime("2015-03-20 {}".format((x.split()[1])), "%Y-%m-%d %H:%M:%S.%f")
    time.append(dt)
    time_julian.append((Time(dt)).jd)
    raw_rv.append(float(x.split()[4]))
for x in range(0,len(lines)):
    UTC_time.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f") + timedelta(seconds=float(50))))
f.close()
#reiners model one residuals and intensity 
with open('model1_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model1 = []
    residual_model1 = []
    for row in cf:
        time_residuals_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model1.append(float(row[1]))
with open('model2_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
with open('new_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_new = []
    residual_new = []
    for row in cf:
        time_residuals_new.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new.append(float(row[1]))
with open('flux_data.csv') as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))     

raw_rv = np.array(raw_rv)
raw_rv -= raw_rv[-1]

GRASS_rv = np.array(GRASS_rv)
GRASS_rv -= GRASS_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb) 
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb)
RV_list_cb -= RV_list_cb[-1]

RV_list_cb_new = np.array(RV_list_cb_new)
RV_list_cb_new -= RV_list_cb_new[-1]

residual_slope = ((raw_rv[46] - RV_list_cb[46])-(raw_rv[0] - RV_list_cb[0]))/46
linear_correction = residual_slope*range(0,46) + (raw_rv[0] - RV_list_cb[0])

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 3, label = "GRASS")
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((GRASS_rv - RV_list_no_cb)**2))/len(GRASS_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_no_cb[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_no_cb[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_new, residual_new, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, (GRASS_rv) - RV_list_no_cb, color = 'b', marker = "x", s = 3) 
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_cb)**2))/len(raw_rv - RV_list_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals       
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb_new, color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_cb_new)**2))/len(raw_rv - RV_list_cb_new)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb_new[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb_new[46:-1], color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_cb_new.png")
plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(UTC_time, intensity_list/intensity_list[129], color = 'r', label = "model")
# ax1.scatter(time_flux_model1, flux_model1, color = 'k', label = "Reiners")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity")
# plt.legend()
# plt.savefig("intensity.png")
# plt.show()

# #barycentric correction test: Sun center vs surface
# vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=51.560583 , longi=9.944333, alt=201, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, vb, color = 'g', label = "BC using center")
# plt.scatter(time, bc_list, color = 'r', label = "BC using surface")
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Barycentric Correction")
# plt.legend()
# plt.savefig("barycentric_correction.png")
# plt.show()

# #extiction test: differential vs absolute 
# f=open("Reiners_Data.txt","r")
# lines=f.readlines()[1:]
# airmass = []
# for x in range(0,len(lines)):
#    airmass.append(float((lines[x].split()[2])))
# f.close()
# ext_coefficient = list(np.linspace(0.4, 0.1, num = len(airmass)))
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, airmass, color = 'b', label = "absolute airmass")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Airmass")
# diff_airmass1 = []
# diff_airmass2 = []
# for i in range(0,len(airmass_matrix)):
#     diff_airmass1.append(1 / file[airmass_matrix[i]][()][1][1])
#     diff_airmass2.append(1 / file[airmass_matrix[i]][()][30][10])
# plt.scatter(time, diff_airmass1, color = 'r', label = "airmass for random patch")
# plt.scatter(time, diff_airmass2, color = 'g', label = "airmass for random patch")
# plt.legend()
# plt.savefig("airmass_test")
# plt.show()
