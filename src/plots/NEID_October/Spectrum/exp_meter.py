import os
import h5py
import csv
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def model_flux_exp(file_data, max_intensity):
    intensity_list = file_data["intensity_list"][()]
    intensity_array = [[]]*len(wavelength)
    for i in range(0,len(wavelength)):
        intensity = (file_data[intensity_list[i]][()])
        intensity_array[i].append(intensity/max_intensity[i])
    return intensity_array

def find_closest_index(arr, value):
    return min(range(len(arr)), key=lambda i: abs(arr[i] - value))

max_intensity = [1212197103544.11, 1212197103544.11, 1221568490391.0015, 1221568490391.0015, 1221568490391.0015,  1221568490391.0015, 1223243571572.7612, 1223243571572.7612, 1223243571572.7612, 1223243571572.7612,  1223243571572.7612,  1223243571572.7612, 1233083890048.2249, 1233896233326.561, 1262251409195.169, 1262251409195.169, 1264078174156.69, 1264078174156.69, 1264078174156.69, 1264078174156.69, 1267906167830.0017, 1267906167830.0017]
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

# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=min(exp_meter_wav), vmax=max((exp_meter_wav)))
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

fig = plt.figure()
ax1 = fig.add_subplot()
for wav_ind in range(0, len(exp_meter_wav)):
    wav = exp_meter_wav[wav_ind]
    color = cmap(norm(wav))

    max_value = max(exp_meter_flux[0][wav_ind][2:-5])
    for i in range(0, len(exp_meter_time)):
        if max(exp_meter_flux[i][wav_ind][2:-5]) > max_value:
            max_value = max(exp_meter_flux[i][wav_ind][2:-5])

    for i in range(0, len(exp_meter_time)):
        ax1.plot(exp_meter_time[i][2:-5], (exp_meter_flux[i][wav_ind][2:-5])/max_value, color = color)

ax1.set_xlabel("hour on 10/14")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_ylabel("relative flux") 
plt.savefig("Eclipse_Figures/ExposureMeter/figures/expmeter_full_wavelength.png")
plt.clf()

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

fig = plt.figure()
ax1 = fig.add_subplot()
intensity_exp_meter = []
start = 0
end = 1000
for j in range(1,7):
    file_data = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_{}000_KSSD.jld2").format(j), "r")
    intensity_array = model_flux_exp(file_data, max_intensity)
    for i, value in enumerate(wavelength):
            color = cmap(norm(value))
            ax1.plot(exp_meter_time_csv[start:end], intensity_array[0][i], color=color)
    start += 1000
    end = int("{}000".format(j+1))
ax1.set_xlabel("hour on 10/14")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_ylabel("relative intensity") 
plt.savefig("Eclipse_Figures/ExposureMeter/figures/model.png")
plt.clf()

# ----------------------------------------
# comparison -  GRASS wavelength 

fig = plt.figure()
ax1 = fig.add_subplot()
intensity_exp_meter = []
start = 0
end = 1000
for j in range(1,7):
    file_data = h5py.File(("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_{}000_KSSD.jld2").format(j), "r")
    intensity_array = model_flux_exp(file_data, max_intensity)
    for i, value in enumerate(wavelength):
            color = cmap(norm(value))
            ax1.plot(exp_meter_time_csv[start:end], intensity_array[0][i], color=color)
    start += 1000
    end = int("{}000".format(j+1))

for ind, value in enumerate(wavelength):
    color = cmap(norm(value))
    wav_ind = find_closest_index(exp_meter_wav, value)
    wav = exp_meter_wav[wav_ind]

    max_value = max(exp_meter_flux[0][wav_ind][2:-5])
    for i in range(0, len(exp_meter_time)):
        if max(exp_meter_flux[i][wav_ind][2:-5]) > max_value:
            max_value = max(exp_meter_flux[i][wav_ind][2:-5])

    for i in range(0, len(exp_meter_time)):
        ax1.plot(exp_meter_time[i][2:-5], (exp_meter_flux[i][wav_ind][2:-5])/max_value, color = color)

ax1.set_xlabel("hour on 10/14")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_ylabel("relative flux") 
plt.savefig("Eclipse_Figures/ExposureMeter/figures/comp.png")
plt.clf()