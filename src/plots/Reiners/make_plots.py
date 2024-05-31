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

GRASS_rv = [78.62939490140276, 81.24843570144355, 81.78939116366995, 83.81297660923343, 87.22614810202072, 82.2296294875392, 84.33657849161013, 85.85621872108679, 90.07289919780098, 91.68677598788669, 87.23657309063464, 89.60862188531381, 92.93329464804121, 94.85318639565293, 92.49205096793102, 95.85508378756256, 94.07736187970148, 96.47534782255313, 97.63117838049142, 102.5113947876279, 103.38957082404526, 104.08016866940086, 104.85028437972504, 103.53865869097446, 108.3444198380847, 109.48610135207684, 111.24837837559856, 111.32904960771137, 111.78525562570555, 116.67954179564754, 115.36401142167712, 119.23858839520241, 116.46292354955762, 119.37618241735724, 120.84640250545806, 122.06389932360769, 125.77522855290341, 123.61950851568734, 124.53714888435292, 126.15351316633296, 130.4565200459777, 133.05371429484194, 133.30278857072037, 138.62628750752796, 135.6910465629485, 137.4640917421016, 139.77834399494856, 137.22621224875587, 137.3203924146658, 128.6537570781281, 110.723047743832, 95.08629641280162, 77.86817964969995, 49.604130888706315, 27.164124661103052, 2.438352817896218, -23.633650951385004, -52.38622716223436, -80.56132071158169, -109.70838706264375, -137.90549248176063, -169.94481847142788, -201.92768115972157, -231.44169840815445, -264.8715920208266, -293.6899448608254, -319.98592766603207, -352.81469132808616, -372.4921760601371, -401.5593895027766, -426.64482212664893, -447.52212932698524, -458.53285073782354, -469.5508533972772, -477.5495168072965, -467.9927257599106, -451.882999797124, -423.53969864786694, -379.6232482116226, -321.00048144351007, -250.9137382555622, -156.83860836477479, -49.7110247413881, 82.25837235921335, 219.72938356203326, 351.61774914456083, 486.50007972607125, 607.7390363571029, 704.976577102915, 786.0586313399419, 843.4297232293169, 893.1155056599839, 926.993756309998, 947.7455837820713, 951.1875841224671, 948.9462346675888, 940.599328555862, 924.1555373687429, 903.4878601907474, 873.6431814979318, 845.6141680893993, 816.7869009707568, 789.5912600845526, 760.7687251746553, 727.3066852459317, 697.0443606403621, 665.7982821945313, 637.2934786440393, 607.828496777189, 574.3590674785478, 547.394299476604, 517.1965582013878, 489.8389425759175, 463.4181534382246, 433.8395104279184, 411.2443401397687, 383.25526329340016, 368.80570177293254, 347.3419391732482, 333.4765042933821, 322.61945332351297, 306.3362400856676, 300.28247216083406, 296.14467706830607, 301.14400051535745, 301.74238109591647, 303.1871253499308, 306.0591271641087, 308.1880475585557, 308.4834919412255, 309.88478770028155, 314.65450581208023, 317.9295618426422, 315.25047967313526, 319.6182089306249, 323.358271655229, 327.86386464342485, 329.8726576516273, 326.4316425630681, 329.80775288599256, 331.01549257073566, 335.96068036149154, 336.79819291856643, 340.2452218823678, 341.38863316148735, 348.24796374108945, 343.04726735429904, 345.8882564836176, 348.1442938288266, 355.02873566001836, 367.2767803239427, 371.3610575578098, 369.4206389670937, 371.1142798528155, 374.9183198653015, 379.8071695760134, 378.82364715811707, 384.2231716728009, 385.3949760409923]
GRASS_no_cb = [5.645192702071402, 6.881145953508354, 10.542522618337022, 11.432456128127283, 9.844779056662333, 10.83151339830789, 11.723989348980519, 17.306181747868415, 16.129066275252853, 17.197472926103796, 18.685012066858103, 19.75313943565771, 18.077887289889333, 19.27378063772797, 22.606533937126798, 24.116355983219155, 27.729907196821454, 25.16912000841445, 25.664878139766934, 37.235815060791445, 27.87274685746384, 31.092891493428322, 30.841591773132755, 33.94494323124, 34.87640994067639, 32.38674387094685, 40.82943021669355, 40.106063919250325, 44.09494660822781, 41.6377971413114, 48.02602708943018, 45.10785419464569, 48.61074556393827, 48.12174316374793, 53.23646037774674, 54.023170496829266, 53.79848681932704, 56.943522674810936, 56.54168971099661, 58.45321199215348, 58.01202587011057, 61.377554838387745, 62.2228920763322, 62.93516973551186, 65.22636220433539, 69.8539654950621, 68.9294321937821, 74.24205492477462, 65.79612134816851, 60.118524326647716, 34.05450814442577, 17.612221749824577, -3.1835721933862065, -25.556258271976727, -49.85948642114956, -76.36412704929235, -104.75639045992361, -128.79975229078485, -161.71662594853535, -188.7434553546292, -213.86262555756574, -245.32257980048078, -276.7296153303542, -305.6865680616393, -327.81988474103196, -352.88643490561225, -381.7117723538565, -408.7832781770904, -429.5192222267998, -458.16131530640854, -475.64788055497394, -484.9505318516713, -496.1297139280524, -494.35484219779204, -492.23336296208146, -477.7486463877481, -458.58420650786184, -420.43609077306786, -369.3516053002866, -302.0314626031023, -219.62308082318222, -118.72920024539886, 2.3849356093391605, 124.57649482904806, 266.1516310377971, 398.68127975168915, 526.2503755765954, 638.0438376703406, 735.3629879167881, 815.9012053451243, 872.7534577689127, 912.8520627337027, 934.7406338160093, 944.6434189390338, 933.8964392471622, 925.4226083943528, 912.1756612998882, 886.8919395305418, 865.4588043480487, 837.2641363504904, 797.1891888969993, 762.1956134655309, 730.9500762339586, 697.4281504641231, 663.6467875470454, 637.9116173498799, 603.2575508577888, 568.5679630165641, 537.1751271005222, 500.6040177286473, 473.30715470947905, 437.2291900350188, 411.5710964467239, 385.25922565773595, 356.42250655694556, 335.07754619505386, 315.06961689196953, 296.20755455251054, 276.2749255430753, 264.3003634301785, 248.20098135785884, 238.70682834952697, 230.93134057460534, 226.30638573745154, 230.70628059746076, 229.58125000609172, 236.0973811163111, 232.63845812343283, 236.00613637948243, 238.44272466644233, 241.87036817870023, 246.3989491628101, 247.24434878830985, 249.45992606015176, 250.4708866759064, 256.48219812005556, 254.4312860155402, 261.3569391603697, 261.76067673139255, 262.4246401792904, 268.61326140685674, 269.2955809182762, 269.5958719841846, 269.89361075867606, 272.4437731948016, 277.209126303372, 275.64471518960204, 280.9944870406398, 284.26493654044367, 287.3180390291937, 296.33513832620275, 300.8290505720266, 303.25186689240564, 303.58298436404107, 307.9239427824477, 308.21508670335515, 311.42970776042625, 311.46789543908335, 316.27328584590555]

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
# plt.savefig("intensity_try.png")
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
