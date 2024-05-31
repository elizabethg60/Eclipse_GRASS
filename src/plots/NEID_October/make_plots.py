import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#GRASS_rv: RV calculated from line profiles + granulation affects on
#GRASS_no_cb: RV calculated from line profiles + granulation affects off
#RV_list_no_cb: RV calculated form weighted projected velocities + no extra CB
#RV_list_cb: RV calculated form weighted projected velocities + with extra CB

line_rv = [10988.483965447942, 10983.031666492592, 10963.656613401143, 10947.888035855463, 10946.950018640888, 10922.22810964551, 10916.687949665025, 10893.404407361471, 10877.806412869024, 10861.894671721493, 10842.476659382484, 10820.368612313652, 10800.333149842123, 10785.845462516243, 10784.833465416144, 10763.904844081779, 10738.701765560661, 10725.38349049003, 10711.434674905804, 10682.944450032415, 10672.955923255391, 10660.377114743787, 10635.71840293295, 10649.368839828878, 10625.62450000527, 10613.905916412437, 10600.019858850632, 10605.896309980213, 10605.731551815461, 10583.865675576708, 10584.089879093573, 10579.165110773465, 10579.710035951532, 10579.213044720445, 10622.210313417898, 10605.733610565105, 10620.728640513698, 10655.37509255263, 10697.608235644411, 10734.76578717747, 10757.175898688985, 10834.365142377308, 10889.651622813471, 10989.514253272338, 11080.867277207957, 11162.90994617398, 11267.025349465668, 11402.471798049606, 11497.612372285743, 11625.35519269941, 11692.23326497383, 11788.656357116011, 11878.739673320017, 11943.16047601283, 11989.72742191598, 12023.326619288446, 12032.428689252818, 12073.87865547019, 12069.194122356008, 12082.183400763492, 12054.103212460712, 12062.763446326113, 12065.322350417824, 12045.19509435767, 12043.141065217524, 12028.015322513793, 12010.951142912525, 11969.723224138917, 11983.014694377238, 11944.999255137536, 11922.8827057398, 11918.20216708876, 11892.405585450724, 11858.259374391211, 11840.309945853058, 11819.836777795223, 11784.71717208622, 11790.769091617143, 11759.39968255211, 11748.034952766557, 11722.707746246137, 11697.026269883203, 11697.039889509053, 11656.343026276782, 11644.134215397535, 11625.070456285577, 11598.17580479622, 11579.529906422913, 11572.505011550653, 11553.690320422638, 11529.160237584718, 11530.789702567972, 11496.730799606396, 11497.211262746281, 11485.702045307291, 11457.267357482877, 11448.094345689924, 11430.025562247527, 11436.992129245958, 11406.7454872159, 11396.427989858008, 11393.768804287087, 11372.617193734875, 11357.38980860548, 11348.9946047973, 11340.14426335226, 11332.227716534819, 11326.798837467482, 11309.09352761164, 11319.842900554777, 11313.771295379258, 11303.093783164491, 11298.982280467502, 11293.878031111111, 11298.344630686039, 11310.930950869651, 11312.361014872857, 11309.214644532864, 11313.59072661132, 11317.79832259729, 11319.197663463103, 11327.188548980856, 11331.046765249377, 11340.768171385416, 11335.260700273704, 11327.564409626426, 11328.133769674454, 11331.452729519739, 11334.351631684187, 11337.602123547107, 11331.449738763375, 11347.144509521017, 11355.798186114684, 11349.043642639706, 11347.17355707626, 11352.493347618161, 11358.072324653494, 11365.22457109068, 11355.583219238577, 11367.505601411545, 11371.835129804767, 11370.22223359933, 11364.707978625345, 11364.59185083627, 11370.788878064188, 11366.9564616391, 11396.705703740414, 11376.431789748061, 11372.173635067968, 11387.614996726239, 11387.903037619877, 11384.686845436974, 11403.254713060696, 11392.829968517333, 11388.81491738289, 11401.7378838769, 11401.149602499281, 11407.049470320108]

#read in data
#GRASS
grass_data = h5py.File("data/neid_all_lines_rv_regular.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
grass_data_no_cb = h5py.File("data/neid_all_lines_rv_off.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
#model
file = h5py.File("data/model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]
RV_list_cb  = file["RV_list_cb"][()]
# intensity_list = file["intensity_list"][()]
#data 
data = pd.read_csv("data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs)
rv_obs -= rv_obs[-1]

for i in range(18,19):#len(lines)):
    GRASS_rv_array = grass_data[GRASS_rv[i]][()]
    GRASS_rv_array = np.array(GRASS_rv_array + vb)
    GRASS_rv_array -= GRASS_rv_array[-1]

    line_rv = np.array(line_rv + vb)
    line_rv -= line_rv[-1]

    GRASS_no_cb_array = grass_data_no_cb[GRASS_no_cb[i]][()]
    GRASS_no_cb_array = np.array(GRASS_no_cb_array + vb)
    GRASS_no_cb_array -= GRASS_no_cb_array[-1]

    RV_list_no_cb_array = file[RV_list_no_cb[i]][()]
    RV_list_no_cb_array = np.array(RV_list_no_cb_array + vb)
    RV_list_no_cb_array -= RV_list_no_cb_array[-1]

    RV_list_cb_array = file[RV_list_cb[i]][()]
    RV_list_cb_array = np.array(RV_list_cb_array + vb)
    RV_list_cb_array -= RV_list_cb_array[-1]

    #rm curve 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs")
    axs[0].scatter(UTC_time, line_rv, color = 'y', marker = "x", s = 18, label = "NEID line RVs") 
    axs[0].plot(UTC_time, RV_list_no_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - No CB")
    axs[0].plot(UTC_time, GRASS_no_cb_array, color = 'b', linewidth = 2, label = "Line RVs - No CB")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    # rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb_array)**2))/len(rv_obs - RV_list_no_cb_array)),2)
    # axs[0].text(UTC_time[-40], 500, "Model RMS {}".format(rms_model_no_cb))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_no_cb_array)**2))/len(rv_obs - GRASS_no_cb_array)),2)
    # axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    rms_model_no_cb = round(np.sqrt((np.nansum((line_rv - RV_list_no_cb_array)**2))/len(line_rv - RV_list_no_cb_array)),2)
    axs[0].text(UTC_time[-40], 500, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv - GRASS_no_cb_array)**2))/len(line_rv - GRASS_no_cb_array)),2)
    axs[0].text(UTC_time[-40], 400, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    # axs[1].scatter(UTC_time, rv_obs - RV_list_no_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, rv_obs - GRASS_no_cb_array, color = 'b', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, line_rv - RV_list_no_cb_array, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, line_rv - GRASS_no_cb_array, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("rm_and_residuals_{}.png".format(lines_no_cb[i]))
    #plt.show()
    plt.clf()

    #rm curve w granulation 
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 18, label = "NEID RVs") 
    axs[0].scatter(UTC_time, line_rv, color = 'y', marker = "x", s = 18, label = "NEID line RVs")
    axs[0].plot(UTC_time, RV_list_cb_array, color = 'r', linewidth = 2, label = "Weighted RVs - CB")
    axs[0].plot(UTC_time, GRASS_rv_array, color = 'b', linewidth = 2, label = "GRASS")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].set_xlabel("Time (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    # rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb_array)**2))/len(rv_obs - RV_list_cb_array)),2)
    # axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv_array)**2))/len(rv_obs - GRASS_rv_array)),2)
    # axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    rms_model_no_cb = round(np.sqrt((np.nansum((line_rv - RV_list_cb_array)**2))/len(line_rv - RV_list_cb_array)),2)
    axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
    rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv - GRASS_rv_array)**2))/len(line_rv - GRASS_rv_array)),2)
    axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    # axs[1].scatter(UTC_time, rv_obs - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, rv_obs - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    axs[1].scatter(UTC_time, line_rv - RV_list_cb_array, color = 'r', marker = "x", s = 3) 
    axs[1].scatter(UTC_time, line_rv - GRASS_rv_array, color = 'b', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time (UTC)", fontsize=12)
    axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("rm_and_residuals_cb_{}.png".format(lines[i]))
    #plt.show()
    plt.clf()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(UTC_time, file[intensity_list[0]][()]/max(file[intensity_list[0]][()]), label = "Model")  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# #plt.legend()
# #plt.savefig("intensity.png")
# plt.show()