from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotlib.pyplot as plt 
mpl = plt.matplotlib 
import numpy as np
import h5py
import csv

path = "/storage/home/efg5335/Desktop/Eclipse_GRASS/src/plots/Reiners/"
#data 
f=open(path + "Reiners_Data.txt","r")
lines=f.readlines()[1:]
raw_time=[]
raw_rv = []
for x in lines[46:-35]:
    raw_time.append(datetime.strptime("2015-03-20 {}".format((x.split()[1])), "%Y-%m-%d %H:%M:%S.%f"))
    raw_rv.append(float(x.split()[4]))
f.close()
#reiners model one residuals  
with open(path + "model1_residuals.csv") as f:
    cf = csv.reader(f)
    time_residuals_model1 = []
    residual_model1 = []
    for row in cf:
        time_residuals_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model1.append(float(row[1]))
#reiners model two residuals  
with open(path + "model2_residuals.csv") as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
#reiners new model one residuals
with open(path + "new_residuals.csv") as f:
    cf = csv.reader(f)
    time_residuals_new1 = []
    residual_new1 = []
    for row in cf:
        time_residuals_new1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new1.append(float(row[1])) 
with open(path + "full_flux_data.csv") as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))  

chi2_array = []
chi2_rv = []
chi2_intensity = []
for sec in range(0,80):
    #model
    file = h5py.File("model_data_{}.jld2".format(sec), "r")
    RV_list_no_cb = file["RV_list_no_cb"][()] 
    RV_list_cb  = file["RV_list_cb"][()]
    RV_list_cb_new  = file["RV_list_cb_new"][()]
    timestamps = file["timestamps"][()]
    intensity_list = file["intensity"][()]
    model_time = []
    for i in timestamps:
        model_time.append(datetime.strptime(str(i).split("'")[1], "%Y-%m-%dT%H:%M:%S.%f"))

    raw_rv = np.array(raw_rv)
    raw_rv -= raw_rv[-1]

    RV_list_no_cb = np.array(RV_list_no_cb) 
    RV_list_no_cb -= RV_list_no_cb[-1]

    RV_list_cb = np.array(RV_list_cb)
    RV_list_cb -= RV_list_cb[-1]

    RV_list_cb_new = np.array(RV_list_cb_new)
    RV_list_cb_new -= RV_list_cb_new[-1]

    f_obs=np.array(raw_rv)
    stf_obs=np.std(raw_rv)
    f_exp=np.array(RV_list_no_cb)
    rv_chi = np.nansum(((f_obs-f_exp)/stf_obs)**2)
    chi2_rv.append(np.nansum(((f_obs-f_exp)/stf_obs)**2))
    
    f_obs=np.array(flux_model1[47:-33])
    stf_obs=np.std(flux_model1[47:-33])
    f_exp=np.array(intensity_list/intensity_list[-1])
    intensity_chi = np.nansum(((f_obs-f_exp)/stf_obs)**2)
    chi2_intensity.append(np.nansum(((f_obs-f_exp)/stf_obs)**2))
    
    chi2_array.append(rv_chi + intensity_chi)

#     #rm curve
#     fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
#     axs[0].scatter(model_time, raw_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
#     axs[0].plot(model_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
#     axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[0].set_xlabel("Time (UTC)")
#     axs[0].set_ylabel("RV [m/s]")
#     axs[0].legend()
#     #residuals
#     axs[1].scatter(model_time, raw_rv - RV_list_no_cb, color = 'r', marker = "x", s = 1)
#     axs[1].scatter(time_residuals_model1[10:-14], residual_model1[10:-14], color = 'k', marker = "x", s = 1)
#     axs[1].scatter(time_residuals_new1[17:-14], residual_new1[17:-14], color = 'b', marker = "x", s = 1)
#     axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[1].set_xlabel("Time (UTC)")
#     axs[1].set_ylabel("Residuals")
#     #plt.savefig("rm_curve_no_cb/{}.png".format(sec))
#     rms_no_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2)
#     axs[0].text(model_time[-20], -400, "RMS {}".format(rms_no_cb))
#     plt.savefig("No_CB/{}_no_cb.png".format(sec))
#     plt.show()

#     #rm curve
#     fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
#     axs[0].scatter(model_time, raw_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
#     axs[0].plot(model_time, RV_list_cb, color = 'r', label = "Model - CB")
#     axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[0].set_xlabel("Time (UTC)")
#     axs[0].set_ylabel("RV [m/s]")
#     axs[0].legend()
#     #residuals
#     axs[1].scatter(model_time, raw_rv - RV_list_cb, color = 'r', marker = "x", s = 1)
#     axs[1].scatter(time_residuals_model2[21:-19], residual_model2[21:-19], color = 'k', marker = "x", s = 1)
#     axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[1].set_xlabel("Time (UTC)")
#     axs[1].set_ylabel("Residuals")
#     rms_cb = round(np.sqrt((np.nansum((raw_rv - RV_list_cb)**2))/len(raw_rv - RV_list_cb)),2)
#     axs[0].text(model_time[-20], -400, "RMS {}".format(rms_cb))
#     plt.savefig("CB/{}_cb.png".format(sec))
#     plt.show()

#     #rm curve
#     fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
#     axs[0].scatter(model_time, raw_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
#     axs[0].plot(model_time, RV_list_cb_new, color = 'r', label = "Model - CB")
#     axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[0].set_xlabel("Time (UTC)")
#     axs[0].set_ylabel("RV [m/s]")
#     axs[0].legend()
#     #residuals
#     axs[1].scatter(model_time, raw_rv - RV_list_cb_new, color = 'r', marker = "x", s = 1)
#     axs[1].scatter(time_residuals_model2[21:-19], residual_model2[21:-19], color = 'k', marker = "x", s = 1)
#     axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     axs[1].set_xlabel("Time (UTC)")
#     axs[1].set_ylabel("Residuals")
#     rms_cb_new = round(np.sqrt((np.nansum((raw_rv - RV_list_cb_new)**2))/len(raw_rv - RV_list_cb_new)),2)
#     axs[0].text(model_time[-20], -400, "RMS {}".format(rms_cb_new))
#     plt.savefig("New_CB/{}_cb_new.png".format(sec))
#     plt.show()
    
#     #intensity
#     fig = plt.figure()
#     ax1 = fig.add_subplot()
#     ax1.plot(model_time, intensity_list/intensity_list[-1], color = 'r', label = "model")
#     ax1.scatter(model_time, flux_model1[47:-33], color = 'k', label = "Reiners")
#     ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#     ax1.set_xlabel("Time (UTC)")
#     ax1.set_ylabel("Relative Intensity")
#     plt.legend()
#     plt.savefig("Intensity/{}_intensity.png".format(sec))
#     plt.show()

min_chi2 = min(chi2_array)
added_sec = range(0, 80)[chi2_array.index(min_chi2)]
print(added_sec)
plt.scatter(range(0, 80), chi2_array)
plt.scatter(added_sec, min_chi2)
plt.xlabel("+ seconds")
plt.ylabel("chi squared")
plt.savefig("chisquare.png")
plt.show()
plt.close()

min_chi2 = min(chi2_intensity)
added_sec = range(0, 80)[chi2_intensity.index(min_chi2)]
print(added_sec)
plt.scatter(range(0, 80), chi2_intensity)
plt.scatter(added_sec, min_chi2)
plt.xlabel("+ seconds")
plt.ylabel("chi squared - intensity")
plt.savefig("chisquare_intensity.png")
plt.show()
plt.close()

min_chi2 = min(chi2_rv)
added_sec = range(0, 80)[chi2_rv.index(min_chi2)]
print(added_sec)
plt.scatter(range(0, 80), chi2_rv)
plt.scatter(added_sec, min_chi2)
plt.xlabel("+ seconds")
plt.ylabel("chi squared - rv")
plt.savefig("chisquare_rv.png")
plt.show()
