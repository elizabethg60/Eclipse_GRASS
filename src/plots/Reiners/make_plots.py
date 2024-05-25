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

#GRASS_rv: RV calculated from line profiles + granulation affects on
#GRASS_no_cb: RV calculated from line profiles + granulation affects off
#RV_list_no_cb: RV calculated form weighted projected velocities + no extra CB
#RV_list_cb: RV calculated form weighted projected velocities + with extra CB

GRASS_rv = [79.09007814646145, 79.97060244578836, 79.72484712701049, 82.14827053720359, 82.89774203495155, 83.35938269151497, 87.98221702831523, 84.80103835132643, 87.43857441187093, 87.30077193826502, 92.8740141626519, 88.36559775064799, 93.97784051236309, 92.31109689079001, 94.82748907706535, 95.14083198799395, 100.76178353674793, 95.22938644988781, 98.29091403344056, 104.51405996939691, 99.96784541741638, 103.45521980726745, 105.45074544364961, 106.89864531430652, 108.53238628673171, 108.37806261394552, 110.02233952819718, 111.60044753985827, 112.06591833905944, 112.09603115224591, 120.00122013834195, 116.18726878866346, 116.71303634386989, 120.73158115811422, 122.01068909405738, 121.4274781298276, 121.396043109978, 120.96694516683536, 129.8327243142946, 129.2227138248606, 132.99846107943807, 130.7810422479657, 132.05497108080337, 135.76839152313548, 137.42518889433398, 137.34908896298188, 139.75177602652315, 141.66279221461994, 135.3830682765045, 127.40179039738436, 112.70876612762379, 92.89163547756581, 77.86186946716894, 54.83111002351205, 24.984729106605823, 2.4795592151616255, -21.284362696741933, -54.77787723375133, -81.056195619501, -110.80307462683903, -140.04814920659777, -171.65351582931902, -201.5309084972803, -229.90827916222526, -258.18188323651816, -291.06289846198945, -317.60256615746914, -347.8128850666609, -379.890858242147, -396.0577160071573, -428.9369674612532, -443.94057968010117, -459.0619101797992, -467.4799231303797, -464.26778052002885, -469.98248957647, -449.7908285082216, -422.20453626182297, -384.9327746556461, -323.4944908302442, -247.60155184286378, -154.42623797080583, -45.104401793492976, 86.83762476617987, 216.50023034346387, 357.41244809297416, 482.74511503629356, 607.6722388546918, 699.2708285110523, 785.4730656434072, 851.2748286756034, 897.6782614118747, 926.1817555629144, 944.2641224004876, 951.6996394443048, 948.8851563343061, 935.8494464392936, 925.8111804308996, 901.1256037897242, 879.9973604452148, 849.9598731562475, 816.4447741689333, 793.4698059344979, 756.4415138611647, 728.0867377831238, 696.3094990593349, 664.3026870757293, 634.3220889341981, 606.0045103984144, 578.3254975978077, 547.0165298958713, 520.3396587853724, 489.0895009818156, 460.7133269039899, 432.4495268556876, 409.8611081799698, 387.5255037861958, 364.33373218246635, 349.7014457797406, 336.3645443127308, 320.3967960822197, 308.4915252666121, 300.7275969622375, 294.96616889517196, 296.2263350079129, 299.659052237558, 302.80578980834053, 304.1419698468446, 307.1257783512524, 309.61238460015966, 308.96863183589966, 311.89021202482405, 315.315169586921, 319.5060781764986, 323.67694916537926, 320.04010184262086, 323.15789068395713, 327.47448042990567, 329.699233980735, 331.988776503402, 329.2735991884892, 337.66845271315765, 336.24877627951696, 340.6335076483043, 342.8654998376139, 344.09333899995954, 345.18172075348195, 351.53239049354045, 350.44855187596875, 356.1441254989089, 364.6972644536705, 368.74937353601854, 371.68344981471336, 370.65358002919896, 373.93560589545774, 377.4009626980851, 376.6463538927171, 376.73776666221084, 386.7643810383523]
GRASS_no_cb = [222.32969569852668, 223.31391782082028, 224.2766494070014, 225.27647028483787, 226.29054157898594, 227.31889074929236, 228.35127344919033, 229.42111714183818, 230.49533366582727, 231.57197628815075, 232.70299244642274, 233.8028944065566, 235.01814533165918, 236.15620000008002, 237.3685573205507, 238.58568591778558, 239.78221168969384, 240.99341039260293, 242.21291391326756, 243.4474320828185, 244.6906520055668, 246.09998470582818, 247.36939016647176, 248.65829292186112, 249.95653735769758, 251.2902521990718, 252.62954802687742, 253.97808871580128, 255.3503896675661, 256.733138229428, 258.1189875781897, 259.5577639406973, 260.97283518195246, 262.4531110432894, 263.8903186201962, 265.3659060888642, 266.87257259144707, 268.3570390125902, 269.84514594490054, 271.4027926124625, 272.93735084573234, 274.4727859773826, 276.0674933784481, 277.637388835568, 279.21709799233116, 280.7935803742216, 282.40109995087687, 284.0337460181637, 280.8187710284141, 272.40975984739396, 256.2784210521819, 240.02201749597205, 221.2338486196744, 199.7519104876716, 177.2537573740183, 152.63890614816336, 121.88128424125749, 94.20515138805193, 65.60958185657667, 34.853059634796494, 4.7436511665798085, -25.866382478121228, -56.52368900764544, -87.55495399107909, -118.33481744180611, -148.21623059761635, -178.47705133459544, -207.25045484730487, -234.2164048402771, -259.74084818586226, -283.04651912747437, -303.2936278023268, -324.7134359134128, -335.6787119158972, -340.1770492689806, -337.39676129309555, -319.6771748828134, -295.3606658082177, -255.9691600904402, -201.56751441085035, -128.32402087498696, -36.98437523835563, 72.66293050229285, 199.88229056734025, 334.1589549979754, 470.0271439002796, 602.550942027546, 720.8689802923923, 822.6997810367062, 906.3906132066132, 970.9861823287616, 1018.212770450472, 1053.462844257592, 1072.378493474741, 1079.5672476696927, 1079.0394142264151, 1071.1755765877804, 1057.81958853137, 1039.9778872768782, 1018.619334659675, 982.6203427642523, 955.8617409348609, 927.7302246109135, 897.9907169401414, 869.2777191430508, 839.3291082426022, 808.5245371741166, 778.1747348200902, 749.2823995287669, 718.7736024239339, 689.7516321493968, 662.3387541765674, 633.5665408307717, 606.6130013574048, 582.0346733264192, 558.134530971535, 531.8809194153238, 511.3399030405232, 492.7397933712887, 476.1943203413273, 461.9217830308568, 450.6420517715402, 442.67643631861364, 438.86702711373084, 440.2820727669156, 442.50257501970316, 444.7381646621715, 446.9816470677119, 449.27869206995035, 451.5000301909343, 453.73815208746913, 455.9721806804671, 458.2041168058997, 460.4385005054155, 462.6946488315879, 464.939405270037, 467.17316409068764, 469.44683241061165, 471.7051423950605, 473.93849122391924, 476.2169011820965, 478.4546053077448, 480.7316504597737, 482.972101280592, 485.31408032208833, 487.56313292284995, 489.8414756716807, 492.09770818055347, 494.33048219380066, 496.6026924587235, 508.64013623893317, 510.9322446275832, 513.1612301938554, 515.3164302780951, 517.6176419344054, 519.8792550968147, 522.0954419555958, 524.326109591872, 526.5731828703124]

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

GRASS_no_cb = np.array(GRASS_no_cb)
GRASS_no_cb -= GRASS_no_cb[-1]

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
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Weighted RVs - No CB")
axs[0].plot(UTC_time, GRASS_no_cb, color = 'b', linewidth = 2, label = "Line RVs - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_no_cb)**2))/len(corrected_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-55], -400, "Weighted RVs - No CB RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((corrected_rv  - GRASS_no_cb )**2))/len(corrected_rv  - GRASS_no_cb)),2)
axs[0].text(UTC_time[-55], -500, "Line RVs - No CB RMS {}".format(rms_grass_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time, corrected_rv - RV_list_no_cb, color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_new1, residual_new1, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, corrected_rv - GRASS_no_cb, color = 'b', marker = "x", s = 3) 
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Weighted RVs - CB")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb)**2))/len(corrected_rv - RV_list_cb)),2)
axs[0].text(UTC_time[-55], -400, "Weighted RVs - CB RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((corrected_rv - GRASS_rv)**2))/len(corrected_rv - GRASS_rv)),2)
axs[0].text(UTC_time[-55], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals       
axs[1].scatter(UTC_time, corrected_rv - RV_list_cb, color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, corrected_rv - GRASS_rv, color = 'b', marker = "x", s = 3)
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
