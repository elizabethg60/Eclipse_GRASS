import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

GRASS_rv = [-957.6407050387614, -952.23416482121, -951.3767591469299, -957.5228198873564, -956.0500913470426, -956.4334520576242, -950.0235546408566, -950.7776287356118, -950.3273406033297, -948.7720078586276, -952.0022715989958, -953.4446319505525, -948.048558875819, -950.4815875972771, -944.7365044042803, -950.3739699150658, -949.5958105891948, -947.8443765016298, -948.0828147544905, -947.9234315459737, -946.0441413264148, -945.9156609722065, -940.072267988213, -944.1053574155213, -942.2277009671724, -937.3842286039379, -938.8249193082844, -938.9169790916927, -932.8786867552561, -933.400020404843, -933.9631582000947, -934.1371503259635, -932.3089812720171, -930.6324648498428, -930.7369476207504, -928.5194591601108, -929.5285518362675, -930.2435911631677, -929.9497652922053, -927.9252148592116, -921.738089497927, -925.9717679730974, -924.3431983138896, -922.5174240188875, -921.6136897696907, -923.1924306142977, -922.1934402193941, -1092.9586119388687, -1102.6076289110565, -1125.639584125821, -1138.2865737050495, -1162.9225887985178, -1181.5775615999878, -1207.657536558772, -1224.755143977197, -1251.6965525724352, -1269.341921344596, -1295.6722135078987, -1313.8116296496532, -1340.3198341198129, -1356.3138706483355, -1380.7985438990183, -1407.9933446795735, -1430.5254795141334, -1445.071708268718, -1473.4130545540036, -1489.1132771356047, -1511.9837176199205, -1536.098545517824, -1560.0726027494786, -1572.257806466378, -1597.0763113057326, -1617.0155035643118, -1632.876126481834, -1648.2052718565474, -1671.0070573692894, -1674.0818910666726, -1685.6547147651615, -1695.286990124197, -1684.152233571177, -1684.4177949464802, -1675.1729623965132, -1652.142482390169, -1604.4836072540988, -1562.3915285944827, -1506.088740854738, -1426.6676875420198, -1330.1440211164593, -1208.4122569267495, -1085.7888419823876, -953.3011480584493, -674.644570754855, -581.9215182596539, -469.90033433125376, -373.3549269403287, -328.33409132851835, -270.6487544174198, -213.19073510038916, -184.63386905184203, -170.07803104860855, -149.70329701887042, -133.26553042423123, -118.9623894722335, -119.74174727638797, -121.96955716191046, -125.820028477974, -140.06100920247783, -149.27524163449561, -152.92063550076313, -157.950092186689, -171.81350549865672, -186.95768467806116, -206.2725430579626, -213.740678683174, -232.2338643701338, -244.1991505187002, -264.72747210682155, -277.1541048955806, -293.78406914261296, -314.0057395291762, -320.9567539854917, -343.6884942923018, -364.64307817032704, -374.63729881033305, -397.5806674607968, -404.40667007173863, -422.6199392870508, -435.7111294781792, -614.1843496791362, -610.1984070908541, -612.0467518947662, -610.3979679177977, -603.9129178881329, -607.7652364201504, -605.4992864065548, -601.6089711113312, -594.5901379781579, -593.0099809489709, -591.231272406493, -590.7584274302442, -587.7934809442817, -582.8256944658812, -585.8297736500132, -582.4121394511559, -448.4716834730782, -286.99513045156607, -288.2709606057779, -289.33726152295935]

#read in data
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]
RV_list_cb  = file["RV_list_cb"][()]
intensity_list = file["intensity_list"][()]
vel_no_cb = file["vel_no_cb"][()]
vel_cb = file["vel_cb"][()]
#data 
data_time = np.loadtxt("Boulder_Data_bin.txt")[:, 0]
data_rv = np.loadtxt("Boulder_Data_bin.txt")[:, 1] * (1.565*10**(-6))
time_stamps_data = []
time_julian = []
for i in data_time:
    dt = datetime.fromtimestamp(i) + timedelta(hours=4)
    time_stamps_data.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=39.995380 , longi=-105.262390 , alt=165.23, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

data_rv = np.array(data_rv)
data_rv -= data_rv[-1]

GRASS_rv = np.array(GRASS_rv + vb)
GRASS_rv -= GRASS_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb + vb)
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb + vb)
RV_list_cb -= RV_list_cb[-1]

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(time_stamps_data[47:-20], data_rv[47:-20], color = 'k', marker = "x", s = 18, label = "Heterodyne RVs") 
axs[0].plot(time_stamps_data[47:-20], RV_list_no_cb[47:-20], color = 'r', linewidth = 3, label = "Model - No CB")
axs[0].plot(time_stamps_data[47:-20], GRASS_rv[47:-20], color = 'b', linewidth = 3, label = "GRASS")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
rms_model_no_cb = round(np.sqrt((np.nansum((data_rv[47:-20] - RV_list_no_cb[47:-20])**2))/len(data_rv[47:-20] - RV_list_no_cb[47:-20])),2)
axs[0].text(time_stamps_data[47:-20][-40], -400, "Model RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((data_rv[47:-20] - GRASS_rv[47:-20])**2))/len(data_rv[47:-20] - GRASS_rv[47:-20])),2)
axs[0].text(time_stamps_data[47:-20][-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].legend(fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(time_stamps_data[47:-20], (data_rv[47:-20]) - RV_list_no_cb[47:-20], color = 'r', marker = "x", s = 3) 
axs[1].scatter(time_stamps_data[47:-20], (GRASS_rv[47:-20]) - RV_list_no_cb[47:-20], color = 'b', marker = "x", s = 3)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)", fontsize=12)
axs[1].set_ylabel("Residuals", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(time_stamps_data[47:-20], data_rv[47:-20], color = 'k', marker = "x", s = 18, label = "Heterodyne RVs") 
axs[0].plot(time_stamps_data[47:-20], RV_list_cb[47:-20], color = 'r', linewidth = 3, label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)", fontsize=12)
axs[0].set_ylabel("RV [m/s]", fontsize=12)
rms_model_cb = round(np.sqrt((np.nansum((data_rv[47:-20] - RV_list_cb[47:-20])**2))/len(data_rv[47:-20] - RV_list_cb[47:-20])),2)
axs[0].text(time_stamps_data[47:-20][-40], -400, "Model RMS {}".format(rms_model_cb))
axs[0].legend(fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#residuals
axs[1].scatter(time_stamps_data[47:-20], (data_rv[47:-20]) - RV_list_cb[47:-20], color = 'k', marker = "x", s = 3)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)", fontsize=12)
axs[1].set_ylabel("Residuals", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("rm_and_residuals_cb.png")
plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(time_stamps_data, intensity_list)  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# plt.savefig("intensity.png")
# plt.show()