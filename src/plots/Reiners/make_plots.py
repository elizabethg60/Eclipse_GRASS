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

GRASS_rv = [187.09797902793548, 188.35612571853082, 187.08940722965514, 184.7346739161293, 189.97596949076566, 188.74445841120226, 194.4215502109147, 192.41615335854289, 193.5170356704864, 195.14384355747362, 199.6103859358756, 199.2249911310877, 194.95096031344502, 200.9249793451887, 199.8567300705662, 202.79650935012899, 204.62563536029762, 207.55168354155992, 203.2392771166369, 207.70413361316994, 206.4076298067908, 207.11867748693092, 210.08524684628514, 216.32976242357773, 212.63861094086602, 215.28338951566752, 214.68414474744947, 220.8890196939489, 221.0517961027122, 221.32564919476974, 222.6461238025465, 225.75197851310807, 225.0200862202873, 225.52552257986062, 230.32578355337475, 232.14117982335563, 233.95482742803006, 232.9123484543904, 231.91792007782146, 236.39532040838716, 234.60639372721167, 237.12804639042085, 240.20192689819484, 243.10850520133283, 246.0129131921984, 245.2095526194004, 246.29471476903575, 249.14901516942476, 248.11367265280404, 241.96312753288143, 231.17746200440246, 217.084058661002, 201.50449214230312, 184.07715628813762, 162.2799350476571, 136.98079012145604, 112.2671763968464, 85.0701306750559, 63.225905068739706, 34.19983797998543, 3.99928531334254, -26.928589656911637, -58.78401028453882, -91.64125589310852, -122.39569742439289, -148.23953726826952, -178.88539663645082, -210.4487398822546, -236.38611620328896, -267.5809225618705, -291.3015350397185, -313.0490367293163, -331.69932339281524, -346.70533478398903, -347.8007429943809, -351.0760054742799, -345.64332066946565, -324.5651520957387, -289.9975262906817, -233.88743884286873, -155.33906643569546, -72.45678762758054, 44.68283899066947, 183.38197580190425, 320.86084803133684, 457.0020637229166, 596.3016180934723, 718.8392467276151, 810.0147203361428, 890.6303915881966, 949.3579960644994, 989.9650908058841, 1010.5258545305919, 1025.3736459939519, 1041.6501012991534, 1032.0103180819547, 1015.3448455066007, 999.605014279216, 975.6683775538135, 950.3814427362956, 919.8950125327265, 891.0218455509838, 865.7456384947776, 833.7028085095425, 802.0238633935369, 774.7997141394012, 742.8931046935876, 712.7092435907252, 675.5965980267191, 646.6738440720699, 621.2346223113894, 597.71567858403, 567.0142226684634, 544.0277642843013, 523.4236414484819, 501.1622043066379, 480.7341975575922, 460.72496493988865, 446.11514842879797, 431.49066427804695, 422.0167900283922, 410.8154630749741, 403.2957895255293, 404.78767528244333, 405.7933952892306, 408.99338958847557, 410.01770083964743, 410.1015479689682, 415.1191353727946, 420.22882426448433, 421.7060680770464, 421.2780702683306, 420.1452913410648, 424.77798399200583, 427.8597033208189, 430.4564604707812, 429.48273538072215, 436.77152666459347, 438.0553847191158, 436.8244501847397, 448.82956483241827, 445.5383513184607, 446.08486668469567, 444.03664136678873, 451.17648297807017, 454.5857339619367, 459.4793827852635, 457.31847841348934, 463.08149563144036, 462.3464059124997, 475.98537603998943, 476.69076257842204, 477.7297588351377, 482.01716547961183, 483.62951548104326, 483.02137932037306, 489.4319347698625, 489.01822816026123, 492.87542792393805]
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
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_no_cb[0:46], color = 'r', marker = "x", s = 1)
axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_no_cb[46:-1], color = 'r', marker = "x", s = 1)
#axs[1].scatter(time_residuals_new, residual_new, color = 'k', marker = "x", s = 1)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals.png")
plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
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
#axs[1].scatter(time_residuals_new, residual_new, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, (GRASS_rv) - RV_list_no_cb, color = 'b', marker = "x", s = 3) 
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
# axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Model - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_cb)**2))/len(raw_rv - RV_list_cb)),2)
# axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals       
# axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb[0:46], color = 'r', marker = "x", s = 1)
# axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb[46:-1], color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_cb.png")
# plt.show()

#rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time[0:46], (raw_rv[0:46]) - linear_correction, color = 'k', marker = "x", s = 15)
# axs[0].scatter(UTC_time[46:-1], raw_rv[46:-1], color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb_new, color = 'r', label = "Model - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_cb_new)**2))/len(raw_rv - RV_list_cb_new)),2)
# axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals
# axs[1].scatter(UTC_time[0:46], raw_rv[0:46] - linear_correction - RV_list_cb_new[0:46], color = 'r', marker = "x", s = 1)
# axs[1].scatter(UTC_time[46:-1], raw_rv[46:-1] - RV_list_cb_new[46:-1], color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_cb_new.png")
# plt.show()

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
