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

ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
norm = mcolors.Normalize(vmin=np.min(orders), vmax=np.max(orders))
cmap = plt.get_cmap('coolwarm')

plt.hist(df_fair['slope'], label = 'Fair')
plt.hist(df_good['slope'], label = 'Good')
plt.hist(df_perfect['slope'], label = 'Perfect')
plt.scatter(ext_coeff_array, [5000]*len(ext_coeff_array), color='k')
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp.png")
plt.clf()

for i, value in enumerate(orders):
    color = cmap(norm(value))
    df_order = df_perfect.loc[(df_perfect['order_phys'] == orders[i])]
    plt.hist(df_order['slope'], color = color)
plt.scatter(ext_coeff_array, [25]*len(ext_coeff_array), color="k")
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp_orders.png")
plt.clf()

df_date = df_perfect.loc[(df_perfect['date'] == '2023-10-15')] #day after eclipse
plt.hist(df_date['slope'])
plt.scatter(ext_coeff_array, [20]*len(ext_coeff_array), color="k")
plt.legend()
plt.xlabel("Extinction Coefficient")
plt.savefig("ext_coeff_comp_date.png")