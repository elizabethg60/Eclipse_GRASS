import csv
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

df = pd.read_csv('daily_extinction.csv')
df_perfect = df.loc[(df['num_obs_fit'] > 200)]
df_good = df.loc[(df['num_obs_fit'] > 150)]
df_fair = df.loc[(df['num_obs_fit'] > 100)]

df_mine = pd.read_csv('/storage/home/efg5335/work/GRASS/data/NEID_three_extinction.csv')
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
norm = mcolors.Normalize(vmin=np.min(orders), vmax=np.max(orders))
cmap = plt.get_cmap('coolwarm')

plt.hist(df_fair['slope'], label = 'Fair')
plt.hist(df_good['slope'], label = 'Good')
plt.hist(df_perfect['slope'], label = 'Perfect')
plt.scatter(df_mine['Ext3'], [5000]*len(df_mine['Ext1']), color='k')
plt.scatter(df_mine['Ext2'], [7500]*len(df_mine['Ext1']), color='k')
plt.scatter(df_mine['Ext1'], [10000]*len(df_mine['Ext1']), color='k')
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp.png")
plt.clf()

for i, value in enumerate(orders):
    color = cmap(norm(value))
    df_order = df_perfect.loc[(df_perfect['order_phys'] == orders[i])]
    plt.hist(df_order['slope'], color = color)
plt.scatter(df_mine['Ext3'], [25]*len(df_mine['Ext1']), color="k")
plt.scatter(df_mine['Ext2'], [50]*len(df_mine['Ext1']), color="k")
plt.scatter(df_mine['Ext1'], [100]*len(df_mine['Ext1']), color="k")
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp_orders.png")
plt.clf()

df_date = df_perfect.loc[(df_perfect['date'] == '2023-10-15')] #day after eclipse
plt.hist(df_date['slope'])
plt.scatter(df_mine['Ext3'], [20]*len(df_mine['Ext1']), color="k")
plt.scatter(df_mine['Ext2'], [30]*len(df_mine['Ext1']), color="k")
plt.scatter(df_mine['Ext1'], [40]*len(df_mine['Ext1']), color="k")
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp_date.png")
