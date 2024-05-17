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

GRASS_rv = [79.09007814646145, 79.97060244578836, 79.72484712701049, 82.14827053720359, 82.89774203495155, 83.35938269151497, 87.98221702831523, 84.80103835132643, 87.43857441187093, 87.30077193826502, 92.8740141626519, 88.36559775064799, 93.97784051236309, 92.31109689079001, 94.82748907706535, 95.14083198799395, 100.76178353674793, 95.22938644988781, 98.29091403344056, 104.51405996939691, 99.96784541741638, 103.45521980726745, 105.45074544364961, 106.89864531430652, 108.53238628673171, 108.37806261394552, 110.02233952819718, 111.60044753985827, 112.06591833905944, 112.09603115224591, 120.00122013834195, 116.18726878866346, 116.71303634386989, 120.73158115811422, 122.01068909405738, 121.4274781298276, 121.396043109978, 120.96694516683536, 129.8327243142946, 129.2227138248606, 132.99846107943807, 130.7810422479657, 132.05497108080337, 135.76839152313548, 137.42518889433398, 137.34908896298188, 139.75177602652315, 141.66279221461994, 135.3830682765045, 127.40179039738436, 112.70876612762379, 92.89163547756581, 77.86186946716894, 54.83111002351205, 24.984729106605823, 2.4795592151616255, -21.284362696741933, -54.77787723375133, -81.056195619501, -110.80307462683903, -140.04814920659777, -171.65351582931902, -201.5309084972803, -229.90827916222526, -258.18188323651816, -291.06289846198945, -317.60256615746914, -347.8128850666609, -379.890858242147, -396.0577160071573, -428.9369674612532, -443.94057968010117, -459.0619101797992, -467.4799231303797, -464.26778052002885, -469.98248957647, -449.7908285082216, -422.20453626182297, -384.9327746556461, -323.4944908302442, -247.60155184286378, -154.42623797080583, -45.104401793492976, 86.83762476617987, 216.50023034346387, 357.41244809297416, 482.74511503629356, 607.6722388546918, 699.2708285110523, 785.4730656434072, 851.2748286756034, 897.6782614118747, 926.1817555629144, 944.2641224004876, 951.6996394443048, 948.8851563343061, 935.8494464392936, 925.8111804308996, 901.1256037897242, 879.9973604452148, 849.9598731562475, 816.4447741689333, 793.4698059344979, 756.4415138611647, 728.0867377831238, 696.3094990593349, 664.3026870757293, 634.3220889341981, 606.0045103984144, 578.3254975978077, 547.0165298958713, 520.3396587853724, 489.0895009818156, 460.7133269039899, 432.4495268556876, 409.8611081799698, 387.5255037861958, 364.33373218246635, 349.7014457797406, 336.3645443127308, 320.3967960822197, 308.4915252666121, 300.7275969622375, 294.96616889517196, 296.2263350079129, 299.659052237558, 302.80578980834053, 304.1419698468446, 307.1257783512524, 309.61238460015966, 308.96863183589966, 311.89021202482405, 315.315169586921, 319.5060781764986, 323.67694916537926, 320.04010184262086, 323.15789068395713, 327.47448042990567, 329.699233980735, 331.988776503402, 329.2735991884892, 337.66845271315765, 336.24877627951696, 340.6335076483043, 342.8654998376139, 344.09333899995954, 345.18172075348195, 351.53239049354045, 350.44855187596875, 356.1441254989089, 364.6972644536705, 368.74937353601854, 371.68344981471336, 370.65358002919896, 373.93560589545774, 377.4009626980851, 376.6463538927171, 376.73776666221084, 386.7643810383523]

#read in data
#model
file = h5py.File("data/model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
intensity_list = file["intensity_list"][()]

#data 
f=open("data/Reiners_Data.txt","r")
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

#reiners model residuals and intensity 
with open('data/model2_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
with open('data/new_residuals_model1.csv') as f:
    cf = csv.reader(f)
    time_residuals_new1 = []
    residual_new1 = []
    for row in cf:
        time_residuals_new1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new1.append(float(row[1]))
with open('data/flux_data.csv') as f:
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

residual_slope = ((raw_rv[46] - RV_list_cb[46])-(raw_rv[0] - RV_list_cb[0]))/46
linear_correction = residual_slope*range(0,46) + (raw_rv[0] - RV_list_cb[0])

corrected_rv = list((raw_rv[0:46]) - linear_correction) + list(raw_rv[46:,])

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_no_cb)**2))/len(corrected_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((corrected_rv - GRASS_rv)**2))/len(corrected_rv - GRASS_rv)),2)
axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time, corrected_rv - RV_list_no_cb, color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_new1, residual_new1, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, corrected_rv - GRASS_rv, color = 'b', marker = "x", s = 3)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb)**2))/len(corrected_rv - RV_list_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals       
axs[1].scatter(UTC_time, corrected_rv - RV_list_cb, color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_cb.png")
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
