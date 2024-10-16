import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def movingmedian(x, window):
    result = []
    start_ind = round(window / 2)
    sides = np.floor(window/2)
    
    while start_ind < len(x)-6:
        result.append(np.median(x[int(start_ind - sides):int(start_ind+sides)]))
        start_ind = start_ind + 1
    return result

#swept data
data_swept = np.load('2024-04-08_swept-rv.npz')
data_time = [datetime.fromisoformat(i) for i in data_swept['timestamp']][3000:-3000]
rv = [i for i in data_swept['rv']][3000:-3000]

fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [500, 1]})
axs[0].scatter(data_time, rv, color = 'r')
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].tick_params(left = False, labelleft = False) 
print(len(data_time)) #2567
#remove outliers
movingmedian_flux = movingmedian(rv, 11)
#delete real data points which lie in the edges
new_time = np.array(data_time[6:-6])
new_flux = np.array(rv[6:-6])
#determine outliers using sigma 3
factor = (np.abs(new_flux - np.array(movingmedian_flux)))
good_data = np.ma.masked_where(factor > 100, new_flux)
good_time = np.ma.masked_where(factor > 100, new_time)
axs[0].scatter(good_time, good_data, color = 'b')
df = pd.DataFrame()
df["Time"] = good_time
df["Value"] = good_data
print(len(good_time)) #2555
#binning
df.set_index('Time', inplace=True)
# Taking mean values for a frequency of 1 minute
df_group = df.groupby(pd.Grouper(level='Time', freq='1T'))['Value'].agg('mean')   
df_group.dropna(inplace=True)
df_group = df_group.to_frame().reset_index()
print(len(df_group)) #224
axs[0].scatter(df_group["Time"], df_group["Value"], color = 'g')
plt.savefig("swept_rm_curve_bin.png")

UTC_time_swept = []
for i in df_group["Time"]:
    UTC_time_swept.append((i + timedelta(seconds=30)).strftime("%Y-%m-%dT%H:%M:%S"))
print(UTC_time_swept)


# #dither data
# data_dither = np.load('2024-04-08_dither-rv.npz')
# data_time = [datetime.fromisoformat(i) for i in data_dither['timestamp']][17000:-16000]
# rv = [i for i in data_dither['rv']][17000:-16000]

# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [500, 1]})
# axs[0].scatter(data_time, rv, color = 'r')
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].tick_params(left = False, labelleft = False) 
# print(len(data_time)) #9853
# #remove outliers
# movingmedian_flux = movingmedian(rv, 11)
# #delete real data points which lie in the edges
# new_time = np.array(data_time[6:-6])
# new_flux = np.array(rv[6:-6])
# #determine outliers using sigma 3
# factor = (np.abs(new_flux - np.array(movingmedian_flux)))
# good_data = np.ma.masked_where(factor > 20, new_flux)
# good_time = np.ma.masked_where(factor > 20, new_time)
# axs[0].scatter(good_time, good_data, color = 'b')
# df = pd.DataFrame()
# df["Time"] = good_time
# df["Value"] = good_data
# print(len(good_time)) #9841
# #binning
# df.set_index('Time', inplace=True)
# # Taking mean values for a frequency of 1 minute
# df_group = df.groupby(pd.Grouper(level='Time', freq='1T'))['Value'].agg('mean')   
# df_group.dropna(inplace=True)
# df_group = df_group.to_frame().reset_index()
# print(len(df_group)) #161
# axs[0].scatter(df_group["Time"], df_group["Value"], color = 'g')
# plt.savefig("dither_rm_curve_bin.png")

# UTC_time_dither = []
# for i in df_group["Time"]:
#     UTC_time_dither.append((i + timedelta(seconds=30)).strftime("%Y-%m-%dT%H:%M:%S"))
# print(UTC_time_dither)