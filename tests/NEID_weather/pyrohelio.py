import csv
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

df_pyrrheliometer = pd.read_csv("neid_ljpyrohelio_chv0_20231009.tel.txt", sep=" ", names=["Date", "Photocell", "Irradiance"])
df_pyrrheliometer["DateTime"] = [datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f") for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Day"] = [str(datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f").date()) for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Hour"] = pd.to_datetime(df_pyrrheliometer['DateTime']).dt.hour 
df_Oct = df_pyrrheliometer.loc[(df_pyrrheliometer['Day'] == "2023-10-14") & (df_pyrrheliometer['Hour'] >= 15) & (df_pyrrheliometer['Hour'] < 19)]

df_pyrrheliometer = pd.read_csv("neid_ljpyrohelio_chv0_20240408.tel.txt", sep=" ", names=["Date", "Photocell", "Irradiance"])
df_pyrrheliometer["DateTime"] = [datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f") for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Day"] = [str(datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f").date()) for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Hour"] = pd.to_datetime(df_pyrrheliometer['DateTime']).dt.hour 
df_April = df_pyrrheliometer.loc[(df_pyrrheliometer['Day'] == "2024-04-08") & (df_pyrrheliometer['Hour'] >= 17) & (df_pyrrheliometer['Hour'] < 20)]

# Create a figure and subplots
fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(10, 5))

# Plot the first figure
axs[0].plot(np.array(df_Oct["DateTime"]), np.array(df_Oct['Irradiance']))
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time On 10/14/23 (UTC)", fontsize=12)
axs[0].set_ylabel("Irradiance", fontsize=12)
# Plot the second figure
axs[1].plot(np.array(df_April["DateTime"]), np.array(df_April['Irradiance']))
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time On 4/8/24 (UTC)", fontsize=12)
# Adjust the layout
plt.tight_layout()
plt.savefig("weather.pdf")