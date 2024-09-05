import os
import h5py
import csv
import pickle
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def model_flux_exp(file_data):
    intensity_list = file_data["intensity_list"][()]
    intensity_array = [[]]*len(wavelength)
    for i in range(0,len(wavelength)):
        intensity = (file_data[intensity_list[i]][()])
        intensity_array[i].append(intensity)
    return intensity_array

def find_closest_index(arr, value):
    return min(range(len(arr)), key=lambda i: abs(arr[i] - value))

wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]

# ----------------------------------------
# Exposure Meter - full wavelength

# October Eclipse
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")["filename"]
timestamps_october = timestamps_full_october[15:-150]

exp_meter_time = []
exp_meter_wav = []
exp_meter_flux = []
for j in range(0,len(timestamps_october[0:-25])):
    inputSpectrum = fits.open('{}/{}'.format(path_october, np.array(timestamps_october)[j]))
    jd = Time(inputSpectrum[14].header["JDREF"], format='jd')
    dt = jd.to_datetime()
    time = []
    for i in inputSpectrum[14].data[0][0]:
        time.append(dt + timedelta(seconds=i))
    exp_meter_time.append(time)
    exp_meter_wav.append(inputSpectrum[14].data[1][30:-14])
    exp_meter_flux.append(inputSpectrum[14].data[4][30:-14])

exp_meter_wav = [array[0] for array in exp_meter_wav[0]]
#bluest: 4292.473 reddest: 9300.471

# # Normalize the values to the range [0, 1]
# norm = mcolors.Normalize(vmin=min(exp_meter_wav), vmax=max((exp_meter_wav)))
# # Create a colormap from blue to red
# cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for wav_ind in range(0, len(exp_meter_wav)):
#     wav = exp_meter_wav[wav_ind]
#     color = cmap(norm(wav))

#     for i in range(0, len(exp_meter_time)):
#         ax1.plot(exp_meter_time[i][2:-5], (exp_meter_flux[i][wav_ind][2:-5]), color = color)

# ax1.set_xlabel("hour on 10/14")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_ylabel("flux") 
# plt.savefig("Eclipse_Figures/ExposureMeter/figures/expmeter_full_wavelength.png")
# plt.clf()

# df = pd.DataFrame({'Wavelength': pd.Series(exp_meter_wav), 'Array2': pd.Series([item.strftime("%Y-%m-%dT%H:%M:%S") for sublist in exp_meter_time for item in sublist[2:-5]])})
# # Save the DataFrame to a CSV file
# df.to_csv('Eclipse_Figures/ExposureMeter/data/exposure_meter_data.csv', index=False)

# ----------------------------------------
# model -  GRASS wavelength 

exp_meter_data = pd.read_csv("Eclipse_Figures/ExposureMeter/data/exposure_meter_data.csv")
exp_meter_time_str = exp_meter_data["Array2"]
exp_meter_time_csv = []
for i in exp_meter_time_str: 
    exp_meter_time_csv.append(datetime.strptime(i, "%Y-%m-%dT%H:%M:%S"))

# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=np.min(wavelength), vmax=np.max(wavelength))
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

# fig = plt.figure()
# ax1 = fig.add_subplot()
# intensity_exp_meter = []
# start = 0
# end = 1000
# for j in range(1,8):
#     file_data = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_{}000_KSSD.jld2").format(j), "r")
#     intensity_array = model_flux_exp(file_data)
#     for i, value in enumerate(wavelength):
#             color = cmap(norm(value))
#             ax1.plot(exp_meter_time_csv[start:end], intensity_array[0][i], color=color)
#     start += 1000
#     end = int("{}000".format(j+1))
# ax1.set_xlabel("hour on 10/14")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_ylabel("intensity") 
# ax1.axvline(x=exp_meter_time_csv[6200])
# ax1.axvline(x=exp_meter_time_csv[-1])
# plt.savefig("Eclipse_Figures/ExposureMeter/figures/model.png")
# plt.clf()

# ----------------------------------------
# comparison -  GRASS wavelength 

file_data_end = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_7000_KSSD.jld2").format(j), "r")
intensity_array_end = model_flux_exp(file_data_end)

fig = plt.figure()
ax1 = fig.add_subplot()
intensity_exp_meter = []
start = 0
end = 1000
neid_full_intensity = []
neid_normal_intensity = []
neid_dic = {key: [] for key in wavelength}
model_full_intensity = []
for j in range(1,8):
    file_data = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_{}000_KSSD.jld2").format(j), "r")
    intensity_array = model_flux_exp(file_data)
    for i, value in enumerate(wavelength):
        color = cmap(norm(value))
        normalize_mean = np.mean(intensity_array_end[0][i][200:-1])
        # ax1.plot(exp_meter_time_csv[start:end], intensity_array[0][i]/normalize_mean, color=color)
        if i == 1:
            model_full_intensity.append(intensity_array[0][i]/normalize_mean)
    start += 1000
    end = int("{}000".format(j+1)) 

for ind, value in enumerate(wavelength):
    color = cmap(norm(value))
    wav_ind = find_closest_index(exp_meter_wav, value)
    wav = exp_meter_wav[wav_ind]

    inner = []
    for i in range(0, len(exp_meter_time)):
        inner.append((exp_meter_flux[i][wav_ind][2:-5]))
    neid_full_intensity.append(inner)

for ind, value in enumerate(wavelength):
    color = cmap(norm(value))
    wav_ind = find_closest_index(exp_meter_wav, value)
    wav = exp_meter_wav[wav_ind]
    normalize_mean = np.mean(np.array([item for sublist in neid_full_intensity[ind] for item in sublist])[6200:-1])

    for i in range(0, len(exp_meter_time)):
        ax1.plot(exp_meter_time[i][2:-5], (exp_meter_flux[i][wav_ind][2:-5])/normalize_mean, color = color)
        neid_dic[value].append((exp_meter_flux[i][wav_ind][2:-5])/normalize_mean)
        if ind == 1:
            neid_normal_intensity.append((exp_meter_flux[i][wav_ind][2:-5])/normalize_mean)

# flattened_list_neid = np.array([item for sublist in neid_normal_intensity for item in sublist])
# flattened_list_model = np.array([item for sublist in model_full_intensity for item in sublist])
# chi2_stat = np.sum((flattened_list_neid - flattened_list_model) ** 2 / flattened_list_model)
# ax1.text(exp_meter_time_csv[500], .5, "Chi Squared {}".format(round(chi2_stat,2)))
# ax1.set_xlabel("hour on 10/14")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_ylabel("relative flux") 
# plt.savefig("Eclipse_Figures/ExposureMeter/figures/comp.png")
# plt.clf()

# ----------------------------------------
# extinction reduced chi square

# ext_data = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/extinction_chi_near1.jld2"), "r")

# extinction_chi_data = ext_data["extinction_chi_data"][()]

# ext_coeff = np.linspace(0.09, 0.2, 20)
# my_dict = {}

# for i in range(0, len(extinction_chi_data)):
#     data_wav = ext_data[extinction_chi_data[i]][()]
#     wav_dict = {}
#     for wav in range(0, len(data_wav)):
#         data_time = ext_data[data_wav[wav]][()]
#         inner_dict = {key: [] for key in ext_coeff}
#         for t in range(0,len(data_time)):
#             for ind in range(0, len(ext_coeff)):
#                     inner_dict[ext_coeff[ind]].append(ext_data[data_time[t]][()][ind])
#         wav_dict[wavelength[wav]] = inner_dict
#     my_dict[i] = wav_dict
# with open('Eclipse_Figures/ExposureMeter/data/data_near1.pkl', 'wb') as pickle_file:
#     pickle.dump(my_dict, pickle_file)

with open('Eclipse_Figures/ExposureMeter/data/data_near1.pkl', 'rb') as pickle_file:
    loaded_dict = pickle.load(pickle_file)

# fig = plt.figure()
# ax1 = fig.add_subplot()
model_dic = {key: [] for key in np.linspace(0.09, 0.2, 20)}
for ind in np.linspace(0.09, 0.2, 20):
    start = 0
    end = 1000
    wav_dic = {key: [] for key in wavelength}
    for j in range(1,8):
        for i, value in enumerate(wavelength):
            color = cmap(norm(value))
            normalize_mean = np.mean(loaded_dict[6][value][ind][200:-1])
            # ax1.plot(exp_meter_time_csv[start:end], loaded_dict[j-1][value][ind]/normalize_mean, color=color)
            wav_dic[value].append(loaded_dict[j-1][value][ind]/normalize_mean)
        start += 1000
        end = int("{}000".format(j+1)) 
    model_dic[ind].append(wav_dic)
# ax1.set_xlabel("hour on 10/14")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_ylabel("relative flux") 
# plt.savefig("Eclipse_Figures/ExposureMeter/figures/ext_flux_near1.png")
# plt.clf()

min_chi = {key: [] for key in wavelength}
for i, value in enumerate(wavelength):
    color = cmap(norm(value))
    min_chi_value = 10
    flattened_list_neid = np.array([item for sublist in neid_dic[value] for item in sublist])
    for ind in np.linspace(0.09, 0.2, 20):
        flattened_list_model = np.array([item for sublist in model_dic[ind][0][value] for item in sublist])
        chi2_stat = np.sum((flattened_list_neid - flattened_list_model) ** 2 / flattened_list_model)
        if chi2_stat < min_chi_value:
            min_ind = ind
            min_chi_value = chi2_stat
    min_chi[value] = min_ind
    model = np.array([item for sublist in model_dic[min_chi[value]][0][value] for item in sublist])
    ax1.plot(exp_meter_time_csv, model, color = color)

    if i == 1:
        ax1.text(exp_meter_time_csv[500], .5, "Chi Squared {}".format(round(min_chi_value,2)))
ax1.set_xlabel("hour on 10/14")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_ylabel("relative flux") 
plt.savefig("Eclipse_Figures/ExposureMeter/figures/ext_comp.png")
plt.clf()

# Specify the name of the CSV file
filename = 'NEID_extinction_coefficient.csv'

# Open the file for writing
with open(filename, 'w', newline='') as file:
    # Create a CSV writer object
    writer = csv.writer(file)
    
    # Write the header (optional)
    writer.writerow(['Key', 'Value'])
    
    # Write the dictionary items
    for key, value in min_chi.items():
        writer.writerow([key, value])